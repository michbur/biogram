#' Tests for function ed1_gen (generating encodings whose distances from a given encoding are equal to 1).
#' 
#' 
#' @param enc encoding  

# aa1 = list(`1` = c("g", "a", "p", "v", "m", "l", "i"), 
#            `2` = c("k", "h"), 
#            `3` = c("d", "e"), 
#            `4` = c("f", "r", "w", "y", "s", "t", "c", "n", "q"))
# aa2<-list(`1`=c("a","b","c"),`2`=c("d","e"), `3`="f")
# aa3<-list(`1`=c("a","b","c"),`2`=c("e"), `3`="d")
# aa4<-list(`1`=c("a","b","c"),`2`=c("e"), `3`="f")

# testing the number of elements in encodings
test1<-function(enc){
  l<-length(unlist(enc))
  ed1_enc<-ed1_gen(enc)
  res1<-lapply(ed1_enc, function(x) length(unlist(x))==l)
  if( FALSE %in% res1){
    stop(paste0("error in total number of elements of encoding ",which(res1==FALSE)))
  }else{
    print("all is well")
  }
}

#testing the uniqueness of elements in encodings
test2<-function(enc){
  all_aa_sorted<-sort(unlist(enc))
  ed1_enc<-ed1_gen(enc)
  res2<-lapply(ed1_enc, function(x) sort(unlist(x))==all_aa_sorted )
  if(FALSE %in% res2){
    stop(paste0("error: wrong aa in encoding ",which(res2==FALSE)))
  }else{
    print("all is well")
  }
}


