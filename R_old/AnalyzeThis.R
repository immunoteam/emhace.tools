

#load("analyze_playground/test_long")
#load("analyze_playground/test_tesla")
# load("analyze_playground/tesla_long")
# data = mtx



Analyze.binders = function(NetMHCpan_out_long){
  as.data.frame(t(apply(NetMHCpan_out_long, 1, function(x){
    vect = rep(F,6)
    if (as.numeric(x[4]) < 2 ){vect[1] = T}
    if (as.numeric(x[4]) < 0.5 ){vect[2] = T}
    if (as.numeric(x[6]) < 2 ){vect[3] = T}
    if (as.numeric(x[6]) < 0.5 ){vect[4] = T}
    if (as.numeric(x[7]) < 500 ){vect[5] = T}
    if (as.numeric(x[7]) < 50 ){vect[6] = T}
    return(c(x, vect))
  })), stringsAsFactors = F) -> out
  for(i in 3:7) {out[,i] = as.numeric(out[,i])}
  for(i in 8:13) {out[,i] = as.logical(out[,i])}
  colnames(out) = c(colnames(NetMHCpan_out_long),"Rank_EL_WB", "Rank_EL_SB",  "Rank_BA_WB", "Rank_BA_SB", "Aff_nm_WB", "Aff_nm_SB")
  return(out)
}

#######################################################################################################################################################

Analyze.ordered = function(NetMHCpan_out_long, orderBy = 1){
  rm = Analyze.binders(NetMHCpan_out_long)
  collink = colnames(rm)[8:13]
  ordered = as.data.frame(t(sapply(unique(rm$allele), function(x){
    tmp = subset(out, out$allele == x)
    binders = apply(tmp[,8:13], 2, function(y){
      sum(as.logical(y))
    })
  })), stringsAsFactors = F) 
  colnames(ordered) = collink
  ordered = ordered[order(ordered[,collink[orderBy]], decreasing = T),]
  return(ordered)
  }

#######################################################################################################################################################

Analyze.SeqLogo = function(NetMHCpan_out_long, cutoff_type = 5, cutoff_value = 500, pep_length = 9){
  string = paste0(c("subset(NetMHCpan_out_long, NetMHCpan_out_long$`", colnames(NetMHCpan_out_long)[3:7][cutoff_type], "` < ", cutoff_value, " & nchar(NetMHCpan_out_long$peptide) == ",pep_length,")"), collapse = "")
  tmp = eval(parse(text = string))
  peplist = lapply(unique(tmp[,1]), function(x){
    tmp[tmp[,1] %in% x,2]
  })
  if (any(lengths(peplist) < 5)) {stop("It seems you don't have enough sequences. At least 5 sequence is needed to create SeqLogo")} else{
  names(peplist) = unique(tmp[,1])
  ggplot() + 
    geom_logo(peplist)+ 
    facet_wrap(~seq_group, nrow = ceiling(length(peplist)/4), scales='free_x')}
}

#######################################################################################################################################################

Analyze.Heatmap = function(NetMHCpan_out_long, type = 5){
  data = NetMHCpan_out_long[,1:7]
  hm_in = dcast(data, allele ~ peptide, value.var = colnames(data)[3:7][type])
  rownames(hm_in) = hm_in[,1]
  hm_in = hm_in[,-1]
  breaksList = seq(0, 50000, by = 10)
  pheatmap(log10(hm_in), 
           color = colorRampPalette(brewer.pal(n = 7, name = "RdYlBu"), bias = 8)(length(breaksList)), 
           breaks = breaksList,
           show_colnames = !ncol(hm_in) > 20,
           legend = )
}



# pheatmap(log10(hm_in), 
#          show_colnames = !ncol(hm_in) > 20,
#          legend = )
#logos skálán ábrázolni, custom legent feliratokkal, kitalálni melyik beállíás lesz a jó, log10 értéket mutassa, felirat: 50


#######################################################################################################################################################    

Analyze.peptides = function(NetMHCpan_out_long){
  rm = Analyze.binders(NetMHCpan_out_long)
  t(sapply(unique(rm$peptide), function(x){
    tmp = subset(rm, rm$peptide == x)
    colSums(tmp[,8:13])
  })) -> pepstat
  return(pepstat)
}



# TODO:
  #generate X-mers import
  #HELP file

  #Bias érték optim.
