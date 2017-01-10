#' Generates a list of encodings whose distances from a given encoding are equal to 1.
#' 
#' 
#' @param enc encoding for which all encodings in distance equal to 1 are seached for.


change_2_groups<-function(enc,i,j,l){
  enc_new<-enc
  enc_new[[i]]<-enc[[i]][-j]
  enc_new[[l]]<-c(enc[[l]], enc[[i]][j])
  enc_new
}

add_new_group<-function(enc,i,j){
  enc_new<-enc
  enc_new[[i]]<-enc[[i]][-j]
  temp_list<-list(enc[[i]][j])
  enc_new1<-c(enc_new, temp_list)
  enc_new1
}

remove_group<-function(enc,i,j,l){
  enc_new<-enc
  enc_new[[l]]<-c(enc[[l]], enc[[i]][j])
  enc_new1<-enc_new[-i]
  enc_new1
}

# general function generating a list of encodings
ed1_gen<-function(enc){
  
  ed1_encodings<-list(enc)
  
  k<-length(enc) #number of grups in enc
  i<-1
  while(i<=k){
    
    n_i<-length(enc[[i]]) #number of aa in group i
    j<-1
    while(j<=n_i){
      
      l<-1
      while(l<=k){
        
        if(l!=i){
          
          if(n_i>1){
            
            enc_new<-change_2_groups(enc,i,j,l)
            ed1_encodings<-c(ed1_encodings, list(enc_new))
            
            
          }else{
            enc_new<-remove_group(enc,i,j,l)
            ed1_encodings<-c(ed1_encodings, list(enc_new))
          }
        }  
        l=l+1   
        
      }
      #adding new group:
      if(n_i>2){
        enc_new<-add_new_group(enc,i,j)
        ed1_encodings<-c(ed1_encodings, list(enc_new))
      }
      if(n_i==2 && j==1){
        enc_new<-add_new_group(enc,i,j)
        ed1_encodings<-c(ed1_encodings, list(enc_new))
      }
      j=j+1
      
    }
    i=i+1
    
  }
  ed1_encodings[-1]
  
}