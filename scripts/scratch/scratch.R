# scratch

tmp <- readRDS(files[429]) #SKCM grf 

maxp <- tmp[tmp$drugname == "667880",]

max
effect_p <- unlist(lapply(maxp$effect, function(y) unlist(lapply(y, function(x) unique(x)[2]))))
mean(effect_p)
# -> p-value are inflated...