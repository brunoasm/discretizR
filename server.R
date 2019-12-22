
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
options(shiny.maxRequestSize=1000*1024^2) 
library(plyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(mclust)

#loads matrix classes and functions
source('matrix_functions.R')

#####functions

#does linear model with variable1 as response, and size + morph*species as predictors
do_species_lm = function(log_values,variable1,size_var,semaphoront_var){
  
  right_terms = c()
  if (!is.na(semaphoront_var)){right_terms=c(right_terms,paste('species*',semaphoront_var,sep=''))}
  else {right_terms=c(right_terms, 'species')}
  if (!is.na(size_var)){right_terms=c(right_terms,paste('offset(',size_var,')',sep=''))}
  
  if (variable1 == 'body_size'){left_term = 'calcbodysize'}
  else {left_term = variable1}
  
  formula_str=paste(left_term,'~',paste(right_terms,collapse='+'),sep='')
  linear_model = lm(as.formula(formula_str),data=log_values,na.action = 'na.exclude')
  cat(formula_str)
  print(summary(linear_model))
  return(linear_model)
}

#does linear model with variable1 as response, and size + morph as predictors
get_residuals = function(log_values,variable1,size_var,semaphoront_var){
  
  right_terms = c()
  if (!is.na(size_var)){right_terms=c(right_terms,paste('offset(',size_var,')',sep=''))}
  if (!is.na(semaphoront_var)){right_terms=c(right_terms,semaphoront_var)}
  if (!(length(right_terms))){right_terms='1'}
  
  if (variable1 == 'body_size'){left_term = 'calcbodysize'}
  else {left_term = variable1}
  
  formula_str=paste(left_term,'~',paste(right_terms,collapse='+'),sep='')
  linear_model = lm(as.formula(formula_str),data=log_values,na.action = 'na.exclude')
  cat(formula_str)
  print(summary(linear_model))
  
  #remove semaphoront effect if non significant
  if (!is.na(semaphoront_var) & any(grepl(semaphoront_var,right_terms))){
    linear_model = step(linear_model,scope=list(lower=as.formula(paste('.~.-',semaphoront_var,sep=''))))
  }
  
  return(linear_model)
}




shinyServer(function(input, output) {
  
  #this reactive values that are used in several functions
  v = reactiveValues(files_uploaded = F, 
                     body_selected = c(),
                     color_graph = NULL,
                     pvalues = NULL,
                     corrs = NULL,
                     excluded_chars = c(),
                     remove_invariable = TRUE,
                     treat_as_ordered = TRUE,
                     tnt_min=0,
                     tnt_max=2)

  #this stores all results for each character that cannot be saved in the table with character information itself
  #e. g. linear model results
  char_results = reactiveValues()
  
  
  
  
  
  ###################reactives
  #returns measurements table
  measurements = reactive({
    return(read.csv(input$rawtable$datapath))
  })
  
  #returns char info table
  char_info = reactive({
    charinfo = read.csv(input$chartable$datapath,as.is = T)
    charinfo[charinfo == ''] = NA
    return(charinfo)
  })
  
  #calculates all sizes that are not named "body_size"
  calculate_other_sizes = reactive({
    req(v$files_uploaded)
    sizes = as.character(unique(v$char_info$size_variable))
    withProgress({
      for (sizestr in sizes){
        if (is.na(sizestr)){} 
        else if (sizestr == 'body_size' ){}
        else if (grepl(',',sizestr)){
          vars = strsplit(sizestr,split = ',')[[1]]
          sizename = paste(vars,collapse='')
          v$measurements[sizename] = apply(v$measurements[vars],1,function(x){exp(mean(log(x)))})
        }
      }
      incProgress(amount = 1/(length(sizes)))
    },message = 'Calculating size variables')
  })
  
  #calculates body size and adds to table
  calculate_body_size = reactive({
    vars = v$body_selected
    datatable = v$measurements
    isolate({v$measurements$calcbodysize = apply(datatable[vars],1,function(x){exp(mean(log(x)))})})
  })
  
  #calculates pairwise correlations and p-values
  #THIS CODE WAS REPLACED AN OBSERVE STATEMENT
  #test_all_correlations = reactive({
  #  req(v$all_residuals)
  #  
  #  all_residuals = matrix(data = NA, nrow = dim(v$logvalues)[1], ncol = length(v$residuals))
  #  
  #  for (i in 1:length(v$residuals)){
  #    if (is.null(v$residuals[[i]])){
  #      if (v$char_info$primary_variable[i] == 'body_size'){
  #        all_residuals[,i] = v$logvalues[,'calcbodysize']
  #      } else {
  #        all_residuals[,i] = v$logvalues[,v$char_info$primary_variable[i]]
  #      }
  #      
  #    } else {
  #      all_residuals[as.integer(rownames(v$residuals[[i]]$data)),i] = v$residuals[[i]]$data$residuals
  #    }
  #    
  #  }
  #  
  #  return(get_all_spearman_corr(all_residuals,v$logvalues$species))
  #  
  #})
  
  #get name for currently selected character
  charname = reactive({
    req(input$variable_selected)
    idx = as.integer(input$variable_selected)
    ndigits = ceiling(log10(dim(v$char_info)[1]))
    return(paste('char',sprintf(paste('%0',ndigits,'d',sep=''),idx),sep='_'))
  })
  
  #get internal discretized character matrix
  get_disc_matrix = reactive({
    req(names(char_results))

    charnames = names(char_results)
    charmatrix = new_MorphMatrix()
    for (charname in charnames){
      req(charname %in% names(v$log_and_discrete))
      i = which(charname == charnames)
      if (i %in% v$excluded_chars){next}
      else{
        
        if (char_results[[charname]]$any_dimorphism_interaction){
          
          regression_intercept = char_results[[charname]]$lm$coefficients[[1]]
          fitted_means = char_results[[charname]]$mclust$parameters$mean
          aspect_ratios = sprintf('%.4g',exp(fitted_means + regression_intercept))
          
          if (is.na(v$char_info$sizename[i])){
            state_names = list(paste(v$char_info$primary_variable[i],'measurement is about',aspect_ratios))
            if (v$char_info$primary_variable[i] == 'body_size'){
              isolate({names(state_names) = paste('body_size (geometric mean of ',paste(v$body_selected,collapse=', '),')',sep='')})
            }
            else {
              names(state_names) = paste(v$char_info$primary_variable[i])
            }
            
            
          }
          else if (v$char_info$size_variable[i] == 'body_size'){
            state_names = list(paste(v$char_info$primary_variable[i],'measurement is about',aspect_ratios,'times the geometric mean of',paste(v$body_selected,collapse=', ')))
            names(state_names) = paste(v$char_info$primary_variable[i],'in relation to the geometric mean of',paste(v$body_selected,collapse=', '))
          }
          else if (grepl(',',v$char_info$size_variable[i])){
            state_names = list(paste(v$char_info$primary_variable[i],'measurement is about',aspect_ratios,'times the geometric mean of',v$char_info$size_variable[i]))
            names(state_names) = paste(v$char_info$primary_variable[i],'in relation to the geometric mean of',v$char_info$size_variable[i])
          }
          else {
            state_names = list(paste(v$char_info$primary_variable[i],'measurement is about',aspect_ratios,'times',v$char_info$size_variable[i]))
            names(state_names) = paste(v$char_info$primary_variable[i],'in relation to',v$char_info$size_variable[i])
          }
          
          for (morph in levels(v$log_and_discrete[[v$char_info$morph_variable[i]]])){
            DF = v$log_and_discrete[v$log_and_discrete[[v$char_info$morph_variable[i]]] == morph,c('species',charname)]
            names(DF)[1] = 'taxon'
            names(DF)[2] = 'state'
            temp_state_names = state_names
            names(temp_state_names) = paste(morph,names(temp_state_names),sep=', ')
            charmatrix = append_char(charmatrix,new_char_from_df(DF,type = 'discrete',statenames = temp_state_names))
          }
          
        }
        else{
          
          #create an indicator variable to signal if morph information should be included in charatcer name (i. e. if there is a linear effect of morph)
          morph_variable = v$char_info[as.integer(gsub('^.+_','',charname())),'morph_variable']
          if (any(grepl(morph_variable,names(char_results[[charname()]]$lm$coefficients)))){
            include_morph_variable = T
            morph_coeffs = char_results[[charname()]]$lm$coefficients[grep(morph_variable,names(char_results[[charname()]]$lm$coefficients))]
          }
          else {
            include_morph_variable = F
          }
          
          
          regression_intercept = char_results[[charname()]]$lm$coefficients[[1]]
          
          
          #temp_lm = lm(calcbodysize~sex,data=v$logvalues)
          #temp_means = char_results[['char_11']]$mclust$parameters$mean
          fitted_means = char_results[[charname()]]$mclust$parameters$mean
          aspect_ratios = sprintf('%.4g',exp(fitted_means + regression_intercept))
          
          
          
          if (is.na(v$char_info$sizename[i])){
            if (include_morph_variable){
              state_names = list(paste(v$char_info$primary_variable[i],
                                       'measurement is about',
                                       aspect_ratios,
                                       'in',
                                       levels(v$logvalues[[morph_variable]])[1]
              ))
              for (i in 2:length(levels(v$logvalues[[morph_variable]]))){
                this_coeff = morph_coeffs[grep(levels(v$logvalues[[morph_variable]])[i],names(morph_coeffs))]
                state_names[[1]] = paste(state_names[[1]],
                                         'and',
                                         aspect_ratios * exp(this_coeff),
                                         'in',
                                         levels(v$logvalues[[morph_variable]])[i])
              }
            }
            else{state_names = list(paste(v$char_info$primary_variable[i],'measurement is about',aspect_ratios))}
            
            if (v$char_info$primary_variable[i] == 'body_size'){
              isolate({names(state_names) = paste('body_size (geometric mean of ',paste(v$body_selected,collapse=', '),')',sep='')})
            }
            else {
              names(state_names) = paste(v$char_info$primary_variable[i])
            }
          }
          else if (v$char_info$size_variable[i] == 'body_size'){
            if (include_morph_variable){
              state_names = list(paste(v$char_info$primary_variable[i],
                                       'measurement is about',
                                       aspect_ratios,
                                       'times the geometric mean of',
                                       paste(v$body_selected,collapse=', '),
                                       'in',
                                       levels(v$logvalues[[morph_variable]])[1]
              ))
              for (i in 2:length(levels(v$logvalues[[morph_variable]]))){
                this_coeff = morph_coeffs[grep(levels(v$logvalues[[morph_variable]])[i],names(morph_coeffs))]
                state_names[[1]] = paste(state_names[[1]],
                                         'and',
                                         aspect_ratios * exp(this_coeff),
                                         'times the geometric mean of',
                                         paste(v$body_selected,collapse=', '),
                                         'in',
                                         levels(v$logvalues[[morph_variable]])[i])
              }
            }
            else {state_names = list(paste(v$char_info$primary_variable[i],'measurement is about',aspect_ratios,'times the geometric mean of',paste(v$body_selected,collapse=', ')))}
            names(state_names) = paste(v$char_info$primary_variable[i],'in relation to the geometric mean of',paste(v$body_selected,collapse=', '))
          }
          else if (grepl(',',v$char_info$size_variable[i])){
            if (include_morph_variable){
              state_names = list(paste(v$char_info$primary_variable[i],
                                       'measurement is about',
                                       aspect_ratios,
                                       'times the geometric mean of',
                                       v$char_info$size_variable[i],
                                       'in',
                                       levels(v$logvalues[[morph_variable]])[1]
              ))
              for (i in 2:length(levels(v$logvalues[[morph_variable]]))){
                this_coeff = morph_coeffs[grep(levels(v$logvalues[[morph_variable]])[i],names(morph_coeffs))]
                state_names[[1]] = paste(state_names[[1]],
                                         'and',
                                         aspect_ratios * exp(this_coeff),
                                         'times the geometric mean of',
                                         v$char_info$size_variable[i],
                                         'in',
                                         levels(v$logvalues[[morph_variable]])[i])
              }
            }
            else {state_names = list(paste(v$char_info$primary_variable[i],'measurement is about',aspect_ratios,'times the geometric mean of',v$char_info$size_variable[i]))}
            names(state_names) = paste(v$char_info$primary_variable[i],'in relation to the geometric mean of',v$char_info$size_variable[i])
          }
          else {
            if (include_morph_variable){
              state_names = list(paste(v$char_info$primary_variable[i],
                                       'measurement is about',
                                       aspect_ratios,
                                       'times',
                                       v$char_info$size_variable[i],
                                       'in',
                                       levels(v$logvalues[[morph_variable]])[1]
              ))
              for (i in 2:length(levels(v$logvalues[[morph_variable]]))){
                this_coeff = morph_coeffs[grep(levels(v$logvalues[[morph_variable]])[i],names(morph_coeffs))]
                state_names[[1]] = paste(state_names[[1]],
                                         'and',
                                         aspect_ratios * exp(this_coeff),
                                         'times',
                                         v$char_info$size_variable[i],
                                         'in',
                                         levels(v$logvalues[[morph_variable]])[i])
              }
            }
            else {state_names = list(paste(v$char_info$primary_variable[i],'measurement is about',aspect_ratios,'times',v$char_info$size_variable[i]))}
            names(state_names) = paste(v$char_info$primary_variable[i],'in relation to',v$char_info$size_variable[i])
          }
          
          
          
          
          DF = v$log_and_discrete[,c('species',charname)]
          names(DF)[1] = 'taxon'
          names(DF)[2] = 'state'
          charmatrix = append_char(charmatrix,new_char_from_df(DF,type = 'discrete',statenames = state_names))
        }        
      }
      
    }
    
    return(charmatrix)
  })
  
  #get internal continuous character matrix
  get_cont_matrix = reactive({
    req(names(char_results), v$log_and_resid)
    
    charnames = names(char_results)
    charmatrix = new_MorphMatrix()
    for (charname in charnames){
      req(charname %in% names(v$log_and_resid))
      i = which(charname == charnames)
      if (i %in% v$excluded_chars){next}
      else{
        if (is.na(v$char_info$sizename[i])){
          if (v$char_info$primary_variable[i] == 'body_size'){
            isolate({char_description = paste('body_size (geometric mean of ',paste(v$body_selected,collapse=', '),')',sep='')})
          }
          else {
            char_description = paste(v$char_info$primary_variable[i])
          }
          
        }
        else if (v$char_info$size_variable[i] == 'body_size'){
          char_description = paste(v$char_info$primary_variable[i],'in relation to the geometric mean of',paste(v$body_selected,collapse=', '))
        }
        else if (grepl(',',v$char_info$size_variable[i])){
          char_description = paste(v$char_info$primary_variable[i],'in relation to the geometric mean of',v$char_info$size_variable[i])
        }
        else {
          char_description = paste(v$char_info$primary_variable[i],'in relation to',v$char_info$size_variable[i])
        }
        
        
        
        if (char_results[[charname]]$any_dimorphism_interaction){

          for (morph in levels(v$log_and_resid[[v$char_info$morph_variable[i]]])){
            
            DF = v$log_and_resid[v$log_and_resid[[v$char_info$morph_variable[i]]] == morph,c('species',charname)]
            names(DF)[1] = 'taxon'
            names(DF)[2] = 'state'

            charmatrix = append_char(charmatrix,new_char_from_df(DF,type = 'continuous',charname = paste(morph,', ',char_description,sep='')))
          }
          
        }
        else{
          
          DF = v$log_and_resid[,c('species',charname)]
          names(DF)[1] = 'taxon'
          names(DF)[2] = 'state'
          
          charmatrix = append_char(charmatrix,new_char_from_df(DF,type = 'continuous',charname = char_description))
          
        }        
      }
      
    }
    a=1
    b=1
    return(charmatrix)
    
  })
  
  #get nexus matrix
  get_nexus_text = reactive({
    make_nexus_file(get_disc_matrix(),
                    remove_invariable = v$remove_invariable,
                    treat_as_ordered = v$treat_as_ordered)
    })
  
  #get tnt matrix
  get_tnt_text = reactive({
    
    
    char_matrix = switch(input$tnt_as_continuous,
                         'cont' = get_cont_matrix(),
                         'disc' = get_disc_matrix())
    #debug(make_tnt_file)
    make_tnt_file(char_matrix, 
                  remove_invariable = v$remove_invariable, 
                  treat_as_ordered = v$treat_as_ordered, 
                  rescale_min = v$tnt_min, 
                  rescale_max = v$tnt_max)
  })
  
  
  ####################dynamic UIs
  # button for body size options
  output$body_size_options = renderUI({
    req(v$files_uploaded)
    actionButton('body_size_dialogue','Body size options')
  })
  
  # pop up window for body size options
  observeEvent(input$body_size_dialogue,{
    showModal(
      modalDialog(
        title = 'Choose which variable to use as body size',
        'If more than one chosen, the geometric average will be used',
        fluidRow(),
        actionButton('close_bodysize','DONE'),
        checkboxGroupInput(
          'checkbox_body',
          'Choose variables:',
          choices = names(measurements()),
          selected = v$body_selected),
        easyClose = FALSE,
        footer = helpText('Click DONE above to close this window')
      )
    )
  })
  
  #button that closes pop up
  observeEvent(input$close_bodysize,{
    calculate_body_size()
    removeModal()
  })
  
  # selection box for variables
  output$select_variable = renderUI({
    req(v$files_uploaded)
    char_options = 1:dim(v$char_info)[1]
    names(char_options) = paste(v$char_info$primary_variable,'|',v$char_info$size_variable)
    names(char_options) = gsub(' \\| NA','',names(char_options))
    names(char_options) = sapply(names(char_options),function(x){
      if(nchar(x) > 30) {paste(strtrim(x,30),'...',sep='')}
      else {x}
    })
    names(char_options) = paste(1:length(names(char_options)),names(char_options),sep='-')
    
    selectInput('variable_selected',
                'Choose character to display:',
                choices = char_options
    )
  })
  
  ###CORRELATION TAB UIs
  #selection box for correlated character groups
  output$corr_group_selector = renderUI({
    req(v$corr_groups)
    
    grp_choices = sort(unique(v$corr_groups))
    names(grp_choices) = grp_choices
    names(grp_choices)[which(grp_choices == 0)] = 'Not evaluated'
    
    
    selectInput('corr_group',
                'Choose a group of correlated characters:',
                choices = grp_choices)
  })
  
  #correlation matrix for correlated character groups
  observe({
    req(input$corr_group)
    if(input$corr_group != 0) {
      out_table = v$corrs[which(v$corr_groups == input$corr_group),which(v$corr_groups == input$corr_group)]
      if(length(out_table) > 1){
        out_table = out_table[,1:(dim(out_table)[2]-1)]
        output$corr_table = renderUI({column(6,fluidRow(renderText({'PAIRWISE CORRELATION'})),
                                             fluidRow(renderTable(out_table,
                                                                  rownames = T, 
                                                                  colnames = T,
                                                                  hover = T,
                                                                  bordered = T,
                                                                  align = 'c',
                                                                  na = '    ',
                                                                  digits = 3)))})
      }
      else {output$corr_table = renderUI({})}
    }
    else {output$corr_table = renderUI({})}
  })
  
  #pvalues for correlated character groups
  observe({
    req(input$corr_group)
    if(input$corr_group != 0){
      out_table = v$pvalues[which(v$corr_groups == input$corr_group),which(v$corr_groups == input$corr_group)]
      if(length(out_table) > 1){
        out_table = out_table[,1:(dim(out_table)[2]-1)]
        output$pvalue_table = renderUI({column(6,fluidRow(renderText({'P-VALUES'})),
                                               fluidRow(renderTable(out_table,
                                                                    rownames = T, 
                                                                    colnames = T,
                                                                    hover = T,
                                                                    bordered = T,
                                                                    align = 'c',
                                                                    na = '    ',
                                                                    digits = 3)))})
      }
      else {output$pvalue_table = renderUI({})}
    }
    else {output$pvalue_table = renderUI({})}
  })
  
  #this will show a message in case there are too few degrees of freedom
  observe({
    req(input$corr_group)
    if(input$corr_group == 0) {
      output$corr_table = renderText({'<strong><font color="red">BE CAREFUL. These characters have too few measurements per species. Character correlation cannot be evaluated.</font></strong>'})
    }
  })
  
  #selection panel to choose characters to exclude
  observe({
    req(input$corr_group)
    
    char_idx = which(v$corr_groups == input$corr_group)
    
    names(char_idx) = paste(v$char_info$primary_variable[char_idx],'|',v$char_info$size_variable[char_idx])
    names(char_idx) = gsub(' \\| NA','',names(char_idx))
    names(char_idx) = sapply(names(char_idx),function(x){
      if(nchar(x) > 80) {paste(strtrim(x,80),'...',sep='')}
      else {x}
    })
    names(char_idx) = paste(char_idx,names(char_idx),sep='-')
    
    
    
    output$corr_char_selector = renderUI({column(12,fluidRow(checkboxGroupInput(inputId = 'excluded_chars',
                                                                                label = 'Check any characters you want to exclude from the output:',
                                                                                choices = char_idx,
                                                                                selected = intersect(char_idx, v$excluded_chars),width = '100%')),
                                                 fluidRow(actionButton('excluded_button','EXCLUDE SELECTED FROM OUTPUT')))})
    
    
  })
  
  #OUTPUT UIs
  
  #button to call mrbayes output dialogue
  output$output_for_mrbayes = renderUI({
    req(v$files_uploaded)
    actionButton('mrbayes_dialogue','Download NEXUS')
  })
  
  #dialogue displayed when buttion is pressed
  observeEvent(input$mrbayes_dialogue,{
    showModal(
      modalDialog(
        title = 'Generate NEXUS file',
        fluidPage(
          fluidRow(tags$b('This downloads a NEXUS file compatible with MorphoBank and Mesquite.')),
          fluidRow('Characters marked for exclusion in correlation tab will be removed.'),
          fluidRow('In cases in which there is variation in the degree of dimorphism between species, there will be one character per morph in the output.'),
          fluidRow(checkboxInput('nexus_remove_invariable',label = 'Remove invariable characters',value = v$remove_invariable)),
          fluidRow(checkboxInput('nexus_treat_as_ordered',label = 'Treat character states as ordered',value = v$treat_as_ordered)),
          fluidRow(downloadButton('download_nexus','Download Nexus file'))
        ),
        ##add an upload button
        easyClose = TRUE
      )
    )
  })
  
  #saves options from NEXUS dialogue
  observeEvent(input$nexus_remove_invariable,{
    v$remove_invariable = input$nexus_remove_invariable
  })
  observeEvent(input$nexus_treat_as_ordered,{
    v$treat_as_ordered = input$nexus_treat_as_ordered
  })
  
  #generates nexus output
  output$download_nexus = downloadHandler(filename = 'discretized_matrix.nexus',
                                          content = function(x){write(get_nexus_text(),x)})
  
  #button to call TNT output
  output$output_for_tnt = renderUI({
    req(v$files_uploaded)
    actionButton('tnt_dialogue','Download TNT')
  })
  
  #TNT dialogue
  observeEvent(input$tnt_dialogue,{
    showModal(
      modalDialog(
        title = 'Generate TNT file',
        fluidPage(
          fluidRow(tags$b('This downloads a TNT file.')),
          fluidRow('Characters marked for exclusion in correlation tab will be removed.'),
          fluidRow('In cases in which there is variation in the degree of dimorphism between species, there will be one character per morph in the output.'),
          fluidRow(''),
          fluidRow(column(5,inputPanel(selectInput('tnt_as_continuous',
                                        label = 'Choose data format',
                                        choices = c('Continuous'='cont','Discretized'='disc'),
                                        selected = 'cont',
                                        width = '70%'))),
                   column(7,uiOutput('tnt_secondary_options'))
                   ),
          fluidRow(downloadButton('download_tnt','Download TNT file'))
        ),
        easyClose = TRUE
      )
    )
  })
  
  #TNT options dependent on continuous or discrete characters
  output$tnt_secondary_options = renderUI({
    if (input$tnt_as_continuous == 'disc'){
      fluidPage(
        fluidRow(checkboxInput('tnt_remove_invariable',label = 'Remove invariable characters',value = v$remove_invariable)),
        fluidRow(checkboxInput('tnt_treat_as_ordered',label = 'Treat character states as ordered',value = v$treat_as_ordered))
      )
      
    }
    else if (input$tnt_as_continuous == 'cont'){
      fluidPage(
        fluidRow(tags$b('Select the range to rescale characters:')),
        fluidRow(numericInput('tnt_min','Minimum: ',value=v$tnt_min,min=0,width = '50%')),
        fluidRow(numericInput('tnt_max','Minimum: ',value=v$tnt_max,min=v$tnt_min,width = '50%'))
      )
    }
    
  })
  
  #saves options from TNT dialogue
  observeEvent(input$tnt_remove_invariable,{
    v$remove_invariable = input$tnt_remove_invariable
  })
  observeEvent(input$tnt_treat_as_ordered,{
    v$treat_as_ordered = input$tnt_treat_as_ordered
  })
  observeEvent(input$tnt_min,{
    v$tnt_min = input$tnt_min
  })
  observeEvent(input$tnt_min,{
    v$tnt_max = input$tnt_max
  })
  
  #generates TNT output
  output$download_tnt = downloadHandler(filename = 'matrix.tnt',
                                          content = function(x){write(get_tnt_text(),x)})
  
  
  ##########################observe_events
  
  
  #this saves measurement table to a reactive value
  observeEvent(input$rawtable,{
    v$measurements = measurements()
    v$files_uploaded = v$meas_uploaded & v$char_uploaded
  })
  
  #this saves character info table to a reactive variable
  observeEvent(input$chartable,{
    v$char_info = char_info()
    v$char_info$sizename = gsub(',','',v$char_info$size_variable)
    v$char_info$sizename[v$char_info$sizename == 'body_size'] = 'calcbodysize'
  })
  
  #this records whether both tables were uploaded
  observe({
    v$files_uploaded = !is.null(v$char_info) & !is.null(v$measurements)
  },priority=10)
  
  #this keeps an updated table with log values
  observeEvent(v$measurements,{
    logfy = function(x){
      if (is.double(x)){
        return(log(x))
      } else {
        return(x)
      }
    }
    isolate({
      v$logvalues = as.data.frame(llply(v$measurements, logfy))
      names(v$logvalues) = names(v$measurements)
    })
    
  },priority = 10)
  
  #this saves in a reactive value which variables will calculate body size
  observeEvent(input$checkbox_body,{
    v$body_selected = input$checkbox_body
  })
  
  #calculate all sizes
  observe({
    calculate_other_sizes()
  })
  
  #do all regressions on size, species, dimorphism and species*dimorphism
  observe({
    req(v$files_uploaded,v$logvalues)
    withProgress({
      
      for (i in 1:dim(v$char_info)[1]){
        ndigits = ceiling(log10(dim(v$char_info)[1]))
        charname = paste('char',sprintf(paste('%0',ndigits,'d',sep=''),i),sep='_')
        
        if (is.na(v$char_info$size_variable[i]) | v$char_info$size_variable[i] != 'body_size'){
          req(v$logvalues$calcbodysize)
          a=1
          temp = do_species_lm(v$logvalues,v$char_info$primary_variable[i],v$char_info$sizename[i],v$char_info$morph_variable[i])
          isolate({
            char_results[[charname]]$splm = temp
          })
          
        } else if (v$char_info$size_variable[i] == 'body_size'){
          req(v$logvalues$calcbodysize)
          temp = do_species_lm(v$logvalues,v$char_info$primary_variable[i],'calcbodysize',v$char_info$morph_variable[i])
          isolate({
            char_results[[charname]]$splm = temp
          })
          
        } 
        incProgress(1/dim(v$char_info)[1])
        
      }
      
    },message = 'Doing linear regression on species, sex and size')
    
  })
  
  #test for dimorphism
  observe({
    
    req(length(names(char_results)) > 0)
    
    for (i in 1:length(names(char_results))){
      charname = names(char_results)[i]
      #is any species x sex interaction significant?
      linear_model = char_results[[charname]]$splm
      pvalues = as.data.frame(summary(linear_model)$coefficient)
      names(pvalues)[length(pvalues)] = 'pvalue'
      pvalues$param = rownames(pvalues)
      sig_params = pvalues %>% dplyr::filter(pvalue < 0.05) %>% dplyr::select(param)
      any_dimorphism_interaction = any(grepl(paste(':',v$char_info$semaphoront_var[i],sep=''),sig_params))
      
      char_results[[charname]]$any_dimorphism_interaction = any_dimorphism_interaction
      
    }
    
  }, priority = 2)
  
  #test pairwise correlations of residuals on size, species and sex
  observe({
    
    req(length(names(char_results)) > 0)
    
    charnames = names(char_results)
    
    pvalues = matrix(NA,nrow = length(names(char_results)), ncol = length(names(char_results)))
    corrs = pvalues
    
    pairs = combn(length(names(char_results)):1,2)
    for (j in 1:dim(pairs)[2]){
      char1 = charnames[pairs[,j][1]]
      char2 = charnames[pairs[,j][2]]
      
      tempres = cor.test(resid(char_results[[char1]]$splm), resid(char_results[[char2]]$splm))
      
      pvalues[matrix(pairs[,j],ncol=2)] = tempres$p.value
      corrs[matrix(pairs[,j],ncol=2)] = tempres$estimate
      
      rownames(pvalues) = rownames(corrs) = colnames(pvalues) = colnames(corrs) = 1:dim(corrs)[1]
      
    }
    
    v$pvalues = pvalues
    v$corrs = corrs
  })
  
  #find clusters of highly correlated variables
  observe({
    req(v$corrs, v$pvalues)
    nchars = length(names(char_results))
    grps = rep(-1,nchars)
    
    #start by pulling out variables with too few degrees of freedom for residuals
    #this will be coded as 0
    for (i in 1:nchars){
      charname = names(char_results)[i]
      if (char_results[[charname]]$splm$df.residual < 10){
        grps[i] = 0
      }
    }
    
    #now we will loop through characters again and group non-zero characters that form groups with significant correlation above the threshold
    for (i in 1:nchars){
      charname = names(char_results)[i]
      if (grps[i] == -1){
      
        highly_correlated = c(which(abs(v$corrs[i,]) > input$mincorr),which(abs(v$corrs[,i]) > input$mincorr))
        low_pvalue = c(which(v$pvalues[i,] < 0.05),which(v$pvalues[,i] < 0.05))
        corr_chars = setdiff(intersect(highly_correlated, low_pvalue),which(grps == 0))

        
        
        #if correlated characters are found, do a recursion and find all correlated to each other
        if (length(corr_chars)){
          while(TRUE){
            old_corr_chars = corr_chars
            corr_chars = c()
            for (j in old_corr_chars){
              highly_correlated = c(which(abs(v$corrs[j,]) > input$mincorr),which(abs(v$corrs[,j]) > input$mincorr))
              low_pvalue = c(which(v$pvalues[j,] < 0.05),which(v$pvalues[,j] < 0.05))
              corr_chars = setdiff(c(corr_chars,intersect(highly_correlated, low_pvalue)),which(grps == 0))
            }
            if (setequal(old_corr_chars,corr_chars)) {break}
            else (corr_chars = union(corr_chars, old_corr_chars))
          }
          
        }
        
        grps[corr_chars] = max(1,max(grps)+1)
        
        
      }
    }
    
    #after find clusters, assign singletons to their own cluster
    for (i in 1:nchars){
      if (grps[i] == -1){grps[i] = max(1,max(grps)+1)}
    }
    
    #now, reorder to show larger groups first
    grp_order = unlist(as.list(sort(table(grps),decreasing = T)))
    grp_names = names(grp_order)
    if (0 %in% grps){
      grp_order = 1:(length(grp_order)-1)
      names(grp_order) = grp_names[grp_names != '0']
      grp_order['0'] = 0
    }
    else {
      grp_order = 1:length(grp_order)
      names(grp_order) = grp_names
    }
    
    grps = grp_order[as.character(grps)]
    
    #save grouping results
    v$corr_groups = grps
    for (i in 1:length(grps)){char_results[[names(char_results)[i]]]$corr_group = grps[i]}
  })
  
  #mark characters for exclusion
  observeEvent(input$excluded_button,{
    showNotification('Excluded list updated!',duration = 1, closeButton = FALSE,type='message')
    all_options = which(v$corr_groups == input$corr_group)
    isolate({v$excluded_chars = union(setdiff(v$excluded_chars,all_options),input$excluded_chars)})
  })
  
  #get all residuals for regressions on size and dimorphism
  observe({
    req(v$files_uploaded,v$logvalues,length(names(char_results)) > 0)
    isolate({
      
      resids = as.data.frame(matrix(NA,nrow=dim(v$logvalues)[1],ncol=length(names(char_results))))
      names(resids) = names(char_results)
      
      withProgress({
        for (i in 1:length(names(char_results))) {
          char = v$char_info[i,]
          charname = names(char_results)[i]
          req(is.logical(char_results[[charname]]$any_dimorphism_interaction))
          a=1
          #isolate({
          if (!is.na(char$size_variable) & char$size_variable == 'body_size'){
            req(v$logvalues$calcbodysize)
            if(char_results[[charname]]$any_dimorphism_interaction){
              char_results[[charname]]$lm = get_residuals(v$logvalues,char$primary_variable,'calcbodysize',NA)
            }
            else{
              char_results[[charname]]$lm = get_residuals(v$logvalues,char$primary_variable,'calcbodysize',char$morph_variable)
            }
          } else {
            if(char_results[[charname]]$any_dimorphism_interaction){
              char_results[[charname]]$lm = get_residuals(v$logvalues,char$primary_variable,char$sizename,NA)
            }
            else{
              char_results[[charname]]$lm = get_residuals(v$logvalues,char$primary_variable,char$sizename,char$morph_variable)
            }
            
          }
          #})
          resids[[charname]] = resid(char_results[[charname]]$lm)
          incProgress(1/dim(v$char_info)[1])
        }
      },message = 'Calculating residuals')
    })
    
    v$log_and_resid = cbind(v$logvalues, resids)
    
    
    
  })
  
  #get gaussian mixture model
  observeEvent(v$log_and_resid,{
    
    chars = grep('char_',names(v$log_and_resid),value=T)
    
    #this fits gaussian mixtures
    withProgress(value=0,message = 'fitting Gaussian mixtures',{
      results = llply(chars,function(x){
        residuals = v$log_and_resid[[x]]
        na_res = is.na(residuals)
        incProgress(1/length(chars))
        fitted_model = densityMclust(residuals[!na_res],modelNames='E',G=1:10)
        
        fitted_model$na.omit = na_res
        
        
        return(fitted_model)
        
      })
    })
    
    names(results) = chars
    
    for (i in 1:length(chars)){
      char_results[[chars[i]]]$mclust = results[[i]]
    }
    
    
    idxs = 1:length(chars)
    names(idxs) = chars
    
    
    charstates = as.data.frame(llply(idxs,function(i){
      states = factor(apply(results[[i]]$z,1,function(y){which(y == max(y))}))
      all_states = factor(rep(NA,length(v$log_and_resid[[names(results)[i]]])), levels = levels(states))
      all_states[sort(setdiff(1:length(all_states),which(results[[i]]$na.omit)))] = states
      return(factor(all_states))
    })
    )
    isolate({
      v$log_and_discrete = v$log_and_resid
      v$log_and_discrete[names(charstates)] = charstates
    })
    
    
  })
  
  
  ################text outputs
  observe({
    req(char_results[[charname()]]$any_dimorphism_interaction)
    
    i = as.integer(input$variable_selected)
    
    if (char_results[[charname()]]$any_dimorphism_interaction){
      output$sd_message = renderText({paste('There was an interaction between',v$char_info$morph_variable[i],'and species. The output will have one characcer per',v$char_info$morph_variable[i])})
    }
    else if (is.na(v$char_info$morph_variable[i])){
      output$sd_message = renderText({paste('No morph or sex informed. The output will have only one character for this variable.')})
    }
    else {
      output$sd_message = renderText({paste('No interaction between',v$char_info$morph_variable[i],'and species. The output will have only one character for this variable.')})
    }
    
    req(v$log_and_resid)
    linear_model = get_residuals(v$logvalues,v$char_info$primary_variable[i],v$char_info$sizename[i],v$char_info$morph_variable[i])
    terms = as.character(linear_model$terms)
    left = terms[2]
    right = terms[3:length(terms)]
    message = paste(right, collapse = '+') %>%
      (function(x)paste(left,'~',x,collapse = ' ')) %>%
      (function(x)paste('Fitted linear model: ',x,collapse = ' '))
    
    output$lm_message = renderText({message})
  })
  
  
  
  ###################plots
  #select color variable
  output$select_graph_color = renderUI({
    req(v$files_uploaded)
    selectInput('color_variable',
                'Choose variable to color graphs:',
                choices = c('species','morph','character state')
    )
  })
 
  #plot regression on size
  observe({
    
    
    req(input$variable_selected,v$logvalues,v$char_info$sizename)
    
    i = as.integer(input$variable_selected)
    a=1
    
    if(is.na(v$char_info$sizename[i])){output$regression_size = renderPlot(plot.new())}
    else{
      
      
      color_variable = input$color_variable
      size_label = function(x){if(x == 'calcbodysize'){'body_size'}else{x}}
      
      output$regression_size = renderPlot({
        
        if (color_variable == 'species'){
          ggplot(v$logvalues) +
            geom_point(aes_string(x=v$char_info$sizename[i],y=v$char_info$primary_variable[i], color='species')) +
            xlab(paste('log(',size_label(v$char_info$sizename[i]),')',sep='')) +
            ylab(paste('log(',v$char_info$primary_variable[i],')',sep='')) +
            theme(plot.margin = unit(c(5.5,75.5,5.5,5.5),'pt')) +
            scale_color_discrete(guide = F)
        }
        
        else if (color_variable == 'morph') {
          req(v$char_info$morph_variable[i])
          ggplot(v$logvalues) +
            geom_point(aes_string(x=v$char_info$sizename[i],y=v$char_info$primary_variable[i], color=v$char_info$morph_variable[i])) +
            xlab(paste('log(',size_label(v$char_info$sizename[i]),')',sep='')) +
            ylab(paste('log(',v$char_info$primary_variable[i],')',sep='')) +
            scale_color_discrete(name='Morph    ')
        }
        
        else if (color_variable == 'character state') {
          req(v$log_and_discrete)
          
          ggplot(v$log_and_discrete) +
            geom_point(aes_string(x=v$char_info$sizename[i],
                                  y=v$char_info$primary_variable[i], 
                                  color=charname())) +
            xlab(paste('log(',size_label(v$char_info$sizename[i]),')',sep='')) +
            ylab(paste('log(',v$char_info$primary_variable[i],')',sep='')) +
            scale_color_discrete(name = 'Character\nstate')
        }
      })
    }
  })
  
  #plot residuals per species
  observe({
    
    i = as.integer(input$variable_selected)
    req(v$log_and_resid, input$color_variable)
    color_variable = input$color_variable
    size_label = function(x){if(x == 'calcbodysize'){'body_size'}else{x}}
    
    output$residuals_species = renderPlot({
      
      if (color_variable == 'species'){
        #ndigits = ceiling(log10(dim(v$char_info)[1]))
        #charname = paste('char',sprintf(paste('%0',ndigits,'d',sep=''),i),sep='_')
        
        ggplot(v$log_and_resid) +
          geom_jitter(aes_string(x='species',y=charname(), color='species'),width = 0.15) +
          xlab('Species') +
          ylab(paste('residuals(',v$char_info$primary_variable[i],')',sep='')) +
          theme(plot.margin = unit(c(5.5,75.5,5.5,5.5),'pt')) +
          scale_color_discrete(guide = F) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
      }
      
      else if (color_variable == 'morph') {
        #ndigits = ceiling(log10(dim(v$char_info)[1]))
        #charname = paste('char',sprintf(paste('%0',ndigits,'d',sep=''),i),sep='_')
        
        req(v$char_info$morph_variable[i])
        ggplot(v$log_and_resid) +
          geom_jitter(aes_string(x='species',y=charname(),color=v$char_info$morph_variable[i]),width = 0.15) +
          xlab('Species') +
          ylab(paste('residuals(',v$char_info$primary_variable[i],')',sep='')) +
          scale_color_discrete(name='Morph    ') +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
      }
      
      else if (color_variable == 'character state') {
        req(v$log_and_discrete)
        #ndigits = ceiling(log10(dim(v$char_info)[1]))
        #charname = paste('char',sprintf(paste('%0',ndigits,'d',sep=''),i),sep='_')
        
        temp_data = v$log_and_resid
        temp_data$discrete = v$log_and_discrete[[charname()]]
        
        ggplot(temp_data) +
          geom_jitter(aes_string(x='species',y=charname(),color='discrete'),width = 0.15) +
          xlab('Species') +
          ylab(paste('residuals(',v$char_info$primary_variable[i],')',sep='')) +
          scale_color_discrete(name = 'Character\nstate') +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
      }
    })
  })
  
  #plots the density of Gaussian mixture
  observe({
    
    req(input$variable_selected, v$log_and_resid, v$log_and_discrete)
    
    
    
    
    
    plot_data = data.frame(residuals = v$log_and_resid[[charname()]],
                           discrete = v$log_and_discrete[[charname()]])
    
    
    
    fitted_means = char_results[[charname()]]$mclust$parameters$mean
    names(fitted_means) = 1:length(fitted_means)
    fitted_variance = char_results[[charname()]]$mclust$parameters$variance$sigmasq
    fitted_prob = char_results[[charname()]]$mclust$parameters$pro
    range_residuals = range(plot_data$residuals,na.rm = T)
    
    fitted_df = data.frame(means=fitted_means, 
                           variance=fitted_variance,
                           probability = fitted_prob)
    
    normal_lines = as.data.frame(alply(fitted_df,1,function(x){
      dnorm(x = seq(range_residuals[1],range_residuals[2],length.out = 200),
            mean = x$means, 
            sd = sqrt(x$variance))# * x$prob
    }))
    
    names(normal_lines) =   names(fitted_means)
    
    normal_lines = normal_lines[levels(plot_data$discrete)]
    
    n_points = 200
    normal_lines = data.frame(x=seq(range_residuals[1],range_residuals[2],length.out = n_points),
                              discrete = rep(names(normal_lines),each = n_points),
                              value = as.vector(as.matrix(normal_lines)))
    
    linetype_data = data.frame(x=plot_data$residuals[!is.na(plot_data$residuals)],y=0,lt=rep_len(c('Empirical\ndistribution','Fitted\nGaussian\ncomponent'),sum(!is.na(plot_data$residuals))))
    
    output$densiplot = renderPlot({
      ggplot(plot_data) +
        geom_density(aes(x=residuals, fill=discrete), alpha=0.2) +
        geom_rug(aes(x=residuals, color=discrete)) +
        geom_line(aes(x=x,y=value,color=discrete),linetype='dashed',data=normal_lines) +
        geom_line(aes(x=x,y=y,linetype=lt),alpha=0,data=linetype_data) +
        scale_color_brewer(type='qual',name='Discretized\ncharacter\nstate') +
        scale_fill_brewer(type='qual', guide = F) +
        scale_linetype_manual(values=c('solid','dashed'),guide=guide_legend(title = '',override.aes = list(alpha=1))) +
        xlab('Residuals') +
        ylab('Density') +
        ggtitle('Distribution of discretized characters')
    })
  })
  
  #plots Mclust density
  observe({
    req(input$variable_selected, v$log_and_resid, v$log_and_discrete)
    mclust_model = char_results[[charname()]]$mclust
    output$mclust_density = renderPlot({
      plot(mclust_model,what='density')
    })
  })
  
  
  #plots BIC
  observe({
    req(input$variable_selected, v$log_and_resid, v$log_and_discrete)
    mclust_model = char_results[[charname()]]$mclust
    output$BICplot = renderPlot({
      plot(mclust_model,what='BIC')
    })
  })
  
  #plots qq
  observe({
    req(input$variable_selected, v$log_and_resid, v$log_and_discrete)
    mclust_model = char_results[[charname()]]$mclust
    output$qqplot = renderPlot({
      densityMclust.diagnostic(mclust_model,type='qq')
    })
  })
  
  #plots cdf
  observe({
    req(input$variable_selected, v$log_and_resid, v$log_and_discrete)
    mclust_model = char_results[[charname()]]$mclust
    output$cdfplot = renderPlot({
      densityMclust.diagnostic(mclust_model,type='cdf')
    })
  })
  
  
  
})
