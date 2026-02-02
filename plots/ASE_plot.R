library(ggplot2)
library(cowplot)
library(tidyr);
library(dplyr)
theme_set(theme_cowplot())
drawBoxPlot=function(out,nam,minExp=5)
{
    snp=nam
    tab<-out %>% group_by(SNP,Gene,Condition,CellType,Sample,Allele) %>% summarise(Count=sum(Count)) %>% spread(Allele,Count,fill=0) %>% 
	as.data.frame()
    tab=tab %>% unite(Name,Gene,SNP,sep="_")
    tab["Ratio"]=tab["alt"]/(tab["alt"]+tab["ref"])
    tab["Tot"]=tab["alt"]+tab["ref"]

    tab=tab[tab$Name==nam,]
    alt=strsplit(nam,":")[[1]][4]

    p= ggplot(tab[tab$Tot>minExp,],aes(x=CellType,y=Ratio))+geom_boxplot(outlier.shape = NA)+  #coord_flip()+
	ylab(paste("Proportion phase UMI mapping to alt allele","",sep=""))+geom_hline(yintercept=.5,linetype="dotted")+
	ggtitle(snp)+ylim(c(0,1))+geom_point(position = position_jitter(w = 0.1, h = 0),aes(color=Condition))

    return(p)
}


drawBoxPlotV2 = function(out, nam, region, p.res, minExp = 5) {
    sv = sub("^[^_]*_", "",nam)
    gene = sub("_.*$", "", nam)
    filt.res = p.res %>% filter(Test=='ASE' & Region==region & CellType %in% unique(out$CellType) & Variant == sv & Gene==gene)
    p_vals = setNames(filt.res$padj, filt.res$CellType)

    snp = nam
    tab <- out %>% 
        group_by(SNP, Gene, Condition, CellType, Sample, Allele, CellType_txt) %>% 
        summarise(Count = sum(Count), .groups = "drop") %>% 
        spread(Allele, Count, fill = 0) %>%
        unite(Name, Gene, SNP, sep = "_") %>%
        mutate(
            Ratio = alt / (alt + ref),
            Tot = alt + ref
        ) %>%
        filter(Name == nam, Tot > minExp)
    
    get_stars = function(p) {
        if (is.na(p)) return("")
        if (p <= 0.0001) return("****")
        if (p <= 0.001) return("***")
        if (p <= 0.01) return("**")
        if (p <= 0.05) return("*")
        return("")
    }
    
    anno_df <- data.frame(
        CellType = names(p_vals),
        p_val = as.numeric(p_vals)
    ) %>% mutate(CellType_txt = stringr::str_wrap(CellType, width=14)) %>%
    rowwise() %>%
    mutate(label = get_stars(p_val)) %>%
    filter(CellType_txt %in% unique(tab$CellType_txt)) # Only label types present in plot

    # Plotting
    p = ggplot(tab, aes(x = CellType_txt, y = Ratio)) +
        geom_boxplot(outlier.shape = NA) +
        geom_point(position = position_jitter(w = 0.1, h = 0), aes(color = Condition)) +
        # Add asterisks at the top (y = 1.05)
        geom_text(data = anno_df, aes(x = CellType_txt, y = 1.05, label = label), size = 5) +
	scale_x_discrete(drop=F) +
        ylab("Phase UMI alt allele proportion") +
        geom_hline(yintercept = .5, linetype = "dotted") +
        ggtitle(paste(region, snp)) +
        ylim(c(0, 1.1)) + # Increased limit to fit stars
        scale_color_brewer(palette = "Dark2") +
        theme_bw() +
        theme(axis.title.y = element_text(size = 14, colour = 'black', face='bold'),
	      axis.title.x = element_blank(),
	      axis.text = element_text(size = 12, colour = 'black', face='bold'))

    #return(p)
    ggsave(paste0(nam,".", region,".pdf"), height=4, width=14)
}



region = 'MTG'
out <- qread("ASE.MTG.qs")

out$CellType[out$CellType=='OPCs'] = 'Oligodendrocyte precursors'
out$CellType = as.factor(out$CellType)
ct.keep = c('Astrocytes','GABA neurons','Glu neurons','Microglia','Oligodendrocyte precursors','Oligodendrocytes')
out = out %>% filter(CellType %in% ct.keep)
out = droplevels(out)
levels(out$CellType) = 
	c('Astrocytes','GABAergic neurons','Glutamatergic neurons','Microglia','Oligodendrocyte precursors','Oligodendrocytes')
out$CellType_txt = stringr::str_wrap(out$CellType, width=14)
out$CellType_txt = factor(out$CellType_txt)

res <- readRDS("../integrate_ASE/combined_results.rds")

nams = c("ARL17B_chr17_46237501_DEL_-724","ITGA8_chr10_15531400_INS_60","MAPT-AS1_chr17_45778945_INS_55")
for (nam in nams) {
drawBoxPlotV2(out, nam, region, res)
}


region = 'Midbrain'
out <- qread("ASE.MB.fixed.qs")

out$CellType[out$CellType=='OPCs'] = 'Oligodendrocyte precursors'
out$CellType = as.factor(out$CellType)
ct.keep = c('Astrocytes','GABA neurons','Glu neurons','Microglia','Oligodendrocyte precursors','Oligodendrocytes')
out = out %>% filter(CellType %in% ct.keep)
out = droplevels(out)
levels(out$CellType) =
        c('Astrocytes','GABAergic neurons','Glutamatergic neurons','Microglia','Oligodendrocyte precursors','Oligodendrocytes')
out$CellType_txt = stringr::str_wrap(out$CellType, width=14)
out$CellType_txt = factor(out$CellType_txt)

res <- readRDS("../integrate_ASE/combined_results.rds")

nams = c("MAPT-AS1_chr17_45778945_INS_55")
for (nam in nams) {
drawBoxPlotV2(out, nam, region, res)
}


