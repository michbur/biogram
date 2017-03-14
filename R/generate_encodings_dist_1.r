#' #' Generating encodings with distance equal to 1
#' #' 
#' #' Generates a list of encodings whose distances from a given encoding are equal to 1.
#' #' 
#' #' @param enc encoding for which all encodings in distance equal to 1 are searched for.
#' #' 
#' #' @return a list of encodings.
#' #' 
#' #' @details 
#' #' We assume that the distance between two encodings is equal to 1 if only one element is assigned 
#' #' to different groups in these two encodings.
#' #' 
#' #' The main function to generate the encodings is \code{ed1_gen}. It goes through every single element
#' #' of the list of a given encoding, taking it from the group it originally belonged to and putting it into
#' #' another group (an existing one or a newly created empty one), thus creating a new encoding 
#' #' with the distance 1 to the original one. The functions \code{change_2_groups}, \code{add_new_group} 
#' #' and \code{remove_group} are technical functions which make the changes in the assignments of the elements 
#' #' to the groups. 
#' #' 
#' #' @examples
#' #' \code{aa1} is an encoding in form of a list:
#' #' aa1 = list(`1` = c("g", "a", "p", "v", "m", "l", "i"), 
#' #'         `2` = c("k", "h"),
#' #'         `3` = c("d", "e"), 
#' #'         `4` = c("f", "r", "w", "y", "s", "t", "c", "n", "q"))
#' #' To generate the list of encodings whose distances from \code{aa1} are equal to 1,
#' #' we put \code{aa1} as the argument of the function \code{ed1_gen}:
#' #' aa1_dist1_encodings = ed1_gen(aa1)
#' #' 
#' 
#' 
#' 
#' change_2_groups<-function(enc,i,j,l){
#'   enc_new<-enc
#'   enc_new[[i]]<-enc[[i]][-j]
#'   enc_new[[l]]<-c(enc[[l]], enc[[i]][j])
#'   enc_new
#' }
#' 
#' add_new_group<-function(enc,i,j){
#'   enc_new<-enc
#'   enc_new[[i]]<-enc[[i]][-j]
#'   temp_list<-list(enc[[i]][j])
#'   enc_new1<-c(enc_new, temp_list)
#'   enc_new1
#' }
#' 
#' remove_group<-function(enc,i,j,l){
#'   enc_new<-enc
#'   enc_new[[l]]<-c(enc[[l]], enc[[i]][j])
#'   enc_new1<-enc_new[-i]
#'   enc_new1
#' }
#' 
#' # general function generating a list of encodings
#' ed1_gen<-function(enc){
#'   
#'   ed1_encodings<-list(enc)
#'   
#'   k<-length(enc) #number of grups in enc
#'   i<-1
#'   while(i<=k){
#'     
#'     n_i<-length(enc[[i]]) #number of aa in group i
#'     j<-1
#'     while(j<=n_i){
#'       
#'       l<-1
#'       while(l<=k){
#'         
#'         if(l!=i){
#'           
#'           if(n_i>1){
#'             
#'             enc_new<-change_2_groups(enc,i,j,l)
#'             ed1_encodings<-c(ed1_encodings, list(enc_new))
#'             
#'             
#'           }else{
#'             enc_new<-remove_group(enc,i,j,l)
#'             ed1_encodings<-c(ed1_encodings, list(enc_new))
#'           }
#'         }  
#'         l=l+1   
#'         
#'       }
#'       #adding new group:
#'       if(n_i>2){
#'         enc_new<-add_new_group(enc,i,j)
#'         ed1_encodings<-c(ed1_encodings, list(enc_new))
#'       }
#'       if(n_i==2 && j==1){
#'         enc_new<-add_new_group(enc,i,j)
#'         ed1_encodings<-c(ed1_encodings, list(enc_new))
#'       }
#'       j=j+1
#'       
#'     }
#'     i=i+1
#'     
#'   }
#'   ed1_encodings[-1]
#'   
#' }
#' 
#' get_1ed <- function(enc) {
#'   len_enc <- length(enc)
#'   new_encs <- unlist(unlist(lapply(1L:len_enc, function(single_group_id) {
#'     single_group <- enc[[single_group_id]]
#'     len_group <- length(single_group)
#'     potential_groups <- 1L:(len_enc + 1)
#'     
#'     lapply(1L:len_group, function(single_unigram_id) {
#'       # new encoding with removed single element and added empty group
#'       new_enc <- c(list(single_group[-single_unigram_id]),
#'                    enc[-single_group_id],
#'                    list(c())) 
#'       lapply(potential_groups[-single_group_id], function(destiny_group_id) {
#'         new_enc[[destiny_group_id]] <- c(single_group[single_unigram_id], new_enc[[destiny_group_id]])
#'         tmp <- new_enc[lengths(new_enc) != 0]
#'         # if(is.na(tmp[[1]]))
#'         #   browser()
#'         tmp
#'       })
#'     })
#'   }), recursive = FALSE), recursive = FALSE)
#'   
#'   only_unique <- unique(unlist(lapply(new_encs, function(single_enc) {
#'     pasted_enc <- sort(unlist(lapply(single_enc, function(single_group)
#'       paste0(sort(single_group), collapse = "")
#'     )))
#'     paste0(pasted_enc[order(nchar(pasted_enc))], collapse = "_")
#'   })))
#'   
#'   lapply(strsplit(only_unique, "_"), strsplit, split = "")
#' }
