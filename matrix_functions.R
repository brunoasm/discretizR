
#######character matrix class and its methods
setClass('DiscCharacter',representation(id = 'numeric',states = 'list', taxnames = 'character',statenames='list'))
setClass('ContCharacter',representation(id = 'numeric',states = 'list', taxnames = 'character',charname='character'))

setClass('MorphMatrix',representation(taxnames = 'character',chars='list'))


#makes a new empty matrix
new_MorphMatrix = function(){
  charmat = new('MorphMatrix')
}

#makes a new character from a dataframe containing taxon names and states
new_char_from_df = function(DF, type = c('discrete','continuous'),statenames=NULL, charname=NULL){
  
  if(suppressWarnings(grepl(type,'discrete'))){
    states = factor(as.character(unique(DF$state)))
    
    taxon_states = dlply(DF,.(taxon),function(x){
      y = as.integer(unique(factor(x$state, levels = levels(states))))
      if (length(y) > 1){y = y[!is.na(y)]; return(sort(y))}
      else {return(y)}
    })
    
    taxa = names(taxon_states)
    names(taxon_states) = NULL
    
    thischar = new('DiscCharacter')
    thischar@taxnames = taxa
    thischar@id=as.integer(1)
    attributes(taxon_states) = NULL
    thischar@states=taxon_states
    
    if (is.null(statenames) & !is.null(charname)){
      statenames = levels(states)
      thischar@statenames=list(charname=statenames)
      names(thischar@statenames) = charname
    }
    else if (is.null(statenames) & is.null(charname)){
      statenames = levels(states)
      thischar@statenames=list(charname=statenames)
      names(thischar@statenames) = 'character'
    }
    else {
      thischar@statenames=statenames
    }
    
    
  }
  else if(suppressWarnings(grepl(type,'continuous'))){
    taxon_states = dlply(DF,.(taxon),function(x){
      if(all(is.na(x$state))){
        return(c(NA,NA))
      }
      else{
        return(range(x$state,na.rm = T))
      }
      
    })
    taxa = names(taxon_states)
    names(taxon_states) = NULL
    thischar = new('ContCharacter')
    thischar@taxnames = taxa
    thischar@id=as.integer(1)
    thischar@states=taxon_states
    thischar@charname = charname
    
  }
  else{
    stop('Charatcer should be either continuous or discrete')
  }
  return(thischar)
  
}

#appends a character to the end of a matrix
setGeneric("append_char", function(matobj, ...) {
  standardGeneric("append_char")
})

setMethod("append_char", signature(matobj = "MorphMatrix"), function(matobj, char) {
  morphmatrix = matobj
  
  #set the id of the new characters to the next id available
  max_id = tryCatch({max(sapply(morphmatrix@chars,function(x)x@id))},
                    error = function(e){0})
  
  char@id = max_id + 1
  
  #add any new species, with their state set to NA for the old characters
  new_taxa = setdiff(char@taxnames,morphmatrix@taxnames)
  morphmatrix@taxnames = c(morphmatrix@taxnames, new_taxa)
  if(length(morphmatrix@chars)){
    for(i in 1:length(morphmatrix@chars)){
      if (class(morphmatrix@chars[[i]]) == 'DiscCharacter'){
        morphmatrix@chars[[i]]@states = append(morphmatrix@chars[[i]]@states, as.list(rep(NA, length(new_taxa))))
      }
      else if (class(morphmatrix@chars[[i]]) == 'ContCharacter'){
        morphmatrix@chars[[i]]@states = append(morphmatrix@chars[[i]]@states, as.list(rep(c(NA,NA), length(new_taxa))))
      }
      
    }
  }
  
  #order new character according to taxa in the new matrix
  char@states = char@states[match(morphmatrix@taxnames,char@taxnames)]
  if (class(char) == 'DiscCharacter'){for (i in 1:length(char@states)){if (is.null(char@states[[i]])) {char@states[[i]] = NA}}}
  if (class(char) == 'ContCharacter'){for (i in 1:length(char@states)){if (is.null(char@states[[i]])) {char@states[[i]] = c(NA,NA)}}}
  
  #append character to morphmatrix
  morphmatrix@chars = append(morphmatrix@chars, list(char))
  
  #update character taxnames
  for (i in 1:length(morphmatrix@chars)){
    morphmatrix@chars[[i]]@taxnames = morphmatrix@taxnames
  }
  
  return(morphmatrix)
  
  
  
})


#combines two character matrices
join_matrices = function(mat1,mat2){
  for (char in mat2@chars){
    mat1 = append_char(mat1,char)
  }
  return(mat1)
}

#parses a nexus file into a character matrix object
###################To be done to enable joining new characters to an existing matrix
parse_nexus_file = function(nexus_path){}

#parses a character matrix into nexus-formatted text
make_nexus_file = function(char_matrix, remove_invariable = T, treat_as_ordered=T){
  #first, remove invariable characters, if option checked
  if(remove_invariable){
    keeps = unlist(llply(char_matrix@chars,function(x){length(x@statenames[[1]]) != 1}))
    char_matrix@chars = char_matrix@chars[keeps]
  }
  
  
  
  #make taxa block
  taxa_block = c('BEGIN TAXA;',
                 paste('DIMENSIONS NTAX =',length(char_matrix@taxnames),';'),
                 'TAXLABELS',
                 paste("'",char_matrix@taxnames,"'",collapse = '\n',sep=''),
                 ';\nEND;\n')
  
  #make characters block
  #char labels
  char_labels = paste('[',
                      1:length(char_matrix@chars),
                      ']',
                      paste("'",unlist(llply(char_matrix@chars,function(x)names(x@statenames))),"'",sep=''),
                      sep='')
  
  
  #state labels
  state_labels = paste(1:length(char_matrix@chars),
                       unlist(llply(char_matrix@chars, function(x){paste("'",x@statenames[[1]],"'",sep='', collapse = '\n')})),
                       sep = '\n') %>% 
    paste(collapse=',\n')
  
  
  #matrix
  matrix_text = llply(1:length(char_matrix@taxnames), function(i){
    charline = llply(char_matrix@chars,function(x){
      if (length(x@states[[i]]) > 1){
        paste('(',paste(x@states[[i]]-1,collapse=' '),')',sep='')
      } else {
        if (is.na(x@states[[i]])){'?'}
        else {x@states[[i]]-1}
      }
      
    }) %>% 
      paste(collapse='')
    
    return(unlist(charline))
    
  }) %>%
    unlist()
  
  matrix_w_taxa = paste("'",char_matrix@taxnames,"'",sep='') %>%
    (function(x) paste(x,matrix_text,sep='\t\t\t\t\t'))
  
  uniqchars <- function(x) unique(strsplit(x, "")[[1]]) 
  all_symbols = llply(matrix_text,uniqchars) %>% 
    unlist %>% 
    unique %>% 
    grep(pattern='[0-9]',value = T) %>%
    as.integer %>%
    sort
  
  #assemble character block
  char_block = c('BEGIN CHARACTERS;',
                 paste('DIMENSIONS NCHAR =',length(char_matrix@chars),';'),
                 paste('FORMAT DATATYPE = STANDARD GAP = - MISSING = ? SYMBOLS = "',paste(all_symbols, collapse=' '),'";'),
                 'CHARLABELS',
                 char_labels,
                 ';',
                 'STATELABELS',
                 state_labels,
                 ';',
                 'MATRIX',
                 matrix_w_taxa,
                 ';',
                 'END;\n'
  )
  
  #make assumptions block
  if(treat_as_ordered){
    assump_block = paste('BEGIN ASSUMPTIONS;',
                         paste('TYPESET * UNTITLED = ord: 1 -',length(char_matrix@chars)),
                         ';',
                         'END;',
                         sep = '\n')
  }
  else{
    assump_block = paste('BEGIN ASSUMPTIONS;',
                         paste('TYPESET * UNTITLED = unord: 1 -',length(char_matrix@chars)),
                         ';',
                         'END;',
                         sep = '\n')
  }
  
  
  
  return(c('#NEXUS\n\n',taxa_block,char_block,assump_block))
  
  
}

#parses a charatcer matrix into tnt-formatted text
make_tnt_file = function(char_matrix, remove_invariable = T, treat_as_ordered=T, rescale_min=0, rescale_max=2){
  #first, remove invariable characters, if option checked
  if(remove_invariable){
    keeps = unlist(llply(char_matrix@chars,function(x){
      if(class(x) == 'DiscCharacter'){length(x@statenames[[1]]) != 1}
      else{TRUE}
    }))
    char_matrix@chars = char_matrix@chars[keeps]
  }
  
  #check which characters are continuous
  disc_chars = llply(char_matrix@chars,function(x){
    class(x) == 'DiscCharacter'
  }) %>% 
    unlist %>% 
    which
  
  cont_chars = setdiff(1:length(char_matrix@chars), disc_chars)
  
  #make header
  if (length(cont_chars)){
    tnt_header = paste('nstates cont;',
                       'xread',
                       paste(length(char_matrix@chars),length(char_matrix@taxnames),sep=' '),
                       sep='\n'
    )
  }
  else {
    tnt_header = paste('xread',
                       paste(length(char_matrix@chars),length(char_matrix@taxnames),sep=' '),
                       sep='\n'
    )
  }
  
  
  #make continuous character block
  if(length(cont_chars)){
    #first, pull continuous characters and rescale:
    cont_matrix = new_MorphMatrix()
    cont_matrix@taxnames = char_matrix@taxnames
    for (i in cont_chars){
      thischar = char_matrix@chars[[i]]
      max_st = max(unlist(thischar@states),na.rm=T)
      min_st = min(unlist(thischar@states),na.rm=T)
      rescale = function(x){(x-min_st)/(max_st-min_st)*(rescale_max-rescale_min)+rescale_min}
      thischar@states = llply(thischar@states,rescale)
      cont_matrix = append_char(cont_matrix,thischar)
    }
    
    #now, write matrix
    matrix_text = llply(1:length(cont_matrix@taxnames), function(i){
      charline = llply(cont_matrix@chars,function(x){
        if (all(is.na(x@states[[i]]))){
          return('?')
        } else {
          return(sprintf('%.3f',x@states[[i]]) %>% paste(collapse='-'))
        }
        
      }) %>% 
        paste(collapse=' ')
    })
    
    matrix_text = paste(gsub('[[:punct:][:space:]]','_',cont_matrix@taxnames),
                        matrix_text,
                        sep='\t\t\t')
    
    cont_block = paste(c('&[cont]',
                         matrix_text),
                       collapse='\n')
  }
  else{
    cont_block = ''
  }
  
  #make discrete character block
  if(length(disc_chars)){
    #first, pull discrete characters and rescale:
    disc_matrix = new_MorphMatrix()
    disc_matrix@taxnames = char_matrix@taxnames
    for (i in disc_chars){
      thischar = char_matrix@chars[[i]]
      disc_matrix = append_char(disc_matrix,thischar)
    }
    
    
    matrix_text = llply(1:length(disc_matrix@taxnames), function(i){
      charline = llply(disc_matrix@chars,function(x){
        if (length(x@states[[i]]) > 1){
          paste('[',paste(x@states[[i]]-1,collapse=' '),']',sep='')
        } else {
          if (is.na(x@states[[i]])){'?'}
          else {x@states[[i]]-1}
        }
        
      }) %>% 
        paste(collapse='')
      return(unlist(charline))
    }) %>%
      unlist()
    
    matrix_w_taxa = paste(gsub('[[:punct:][:space:]]','_',char_matrix@taxnames),
                          matrix_text,sep='\t\t\t')
    if (length(cont_chars)){
      disc_block = paste(c('&[num]',
                           matrix_w_taxa),
                         collapse='\n')
    }
    else {
      disc_block = paste(matrix_w_taxa, collapse = '\n')
    }
    
  }
  else{
    disc_block = ''
  }
  
  #make cnames block
  cnames_text = llply(c(cont_chars,disc_chars),function(i){
    thischar = char_matrix@chars[[i]]
    if(class(thischar) == 'DiscCharacter'){
      paste(paste('{',i-1,sep=''),
            gsub('[[:punct:][:space:]]','_',names(thischar@statenames)),
            paste(gsub('\\s','_',thischar@statenames[[1]]),collapse = ' ')
      )
    }
    else {
      paste(paste('{',i-1,sep=''),
            gsub('[[:punct:][:space:]]','_',thischar@charname)
      )
    }
  })
  
  cnames_text = paste(cnames_text,';',sep='')
  cnames_block = paste(c('cnames',
                         cnames_text),
                       sep = '\n')
  
  
  
  #make ccode command (i. e. ordered characters)
  if(treat_as_ordered & length(disc_chars)){
    ccode_block = paste('ccode\t\t\t+\t',
                        paste(length(cont_chars):(length(char_matrix@chars)-1), collapse = ' '),
                        '*;',
                        sep = ' ')
  }
  else {
    ccode_block = ''
  }
  
  #assemble TNT file
  tnt_text = c(tnt_header,
               cont_block,
               disc_block,
               ';',
               cnames_block,
               ';',
               ccode_block,
               '\nproc /;\n;'
  )
  
  return(tnt_text)
  
  
}

