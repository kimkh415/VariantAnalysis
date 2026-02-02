library(GWASTools)
library(Seurat)
library(SNPRelate)
#library(eQTLpipeline)
library(dplyr)
library(scAlleleExpression)
library(purrr)
library(qs)
library(VariantAnnotation)
library(MatrixEQTL)
library(sva)

##Makes GDS file
makeGDS<-function(vcfFile="/stanley/levin_asap_storage/612-eqtl/GenotypeData/Clean_vcf/Combined/comb_new.no.chr.vcf.gz",gdsFile="temp.gds")
{
     return(snpgdsVCF2GDS(vcfFile, gdsFile, verbose=FALSE))
}

##Gets the genomic PCA
getPCA<-function(gdsFile="temp.gds")
{
    gdsobj <- snpgdsOpen(gdsFile)
    snpset <- snpgdsLDpruning(gdsobj,autosome.only=TRUE, maf=0.05, missing.rate=0.05, method="corr", slide.max.bp=10e6,ld.threshold=sqrt(0.1))
    snp.pruned <- unlist(snpset, use.names=FALSE)
    pca <- snpgdsPCA(gdsobj, snp.id=snp.pruned,algorithm="randomized")
    return(pca)
}

##adsf
makeSNPFile<-function(vcf,pseudobulk,snpFile,snpsUse=c())
{
    out=list()
    num=0
    print(num)
    vcffile = open(VcfFile(vcf, yieldSize=100000))
    repeat {
        num=num+1
        print(num)
        vcf = readVcf(vcffile);
        if (nrow(vcf) == 0){break}
        dat=genotypeToSnpMatrix(vcf)
        mat=as(dat$genotype, "numeric")
        mat=mat[colnames(pseudobulk),]
        mn=apply(mat,2,function(x){maf=sum(x)/(2*length(x))})
        mat=mat[,mn>.05 & mn<.95 & !is.na(mn)]
        print(dim(mat))
        out[[num]]=t(mat)
    }
    mat=do.call(rbind,out)
    mat=mat[rownames(mat) %in% snpsUse,]
    system(paste0("echo ",paste(c("snpid",colnames(mat)),collapse="\t")," > ",snpFile))
    write.table(mat,snpFile,col.names=F,sep="\t",quote=F,append=T)
    return(c())
}

##Makes the pseudobulk and saves it to expressFile
makePB<-function(meta,expressFile,genesUse=c(),filterCPM=1,numDS=-1,retFile="")
{
    meta["VCF_Name"]=sub("/stanley/levin_asap_storage/612-eqtl/SingleNuc_data/UpdatedAIData/RunNewData/NewRef/RunPipe_new/samp_","",meta$dirASE)
    
    print("Pseudocbulk meta data")
    meta_PB=meta[!duplicated(meta$VCF_Name),]
    rownames(meta_PB)=meta_PB$VCF_Name

    print("Pseudobulk expression")
    dat=""
    if(!file.exists(retFile))
    {
        dat=LoadExpressionPB(meta)
    }
    else{
        dat=qread(retFile)[[1]]
    }
    if(numDS>0)
    {
        meta_PB=meta_PB[sample(rownames(metaPB),numDS),]
    }
    dat=dat[,rownames(meta_PB)]
    ret=list("data"=dat,"meta"=meta_PB)

    ##normalize
    for(i in colnames(dat)){dat[i]=1000000*dat[,i]/sum(dat[,i])}
    dat=dat[apply(dat,1,function(x) mean(x>filterCPM)>.25),]
    
    print("ComBat normalize")
    datCombat=ComBat(dat, meta_PB[,"batch"], mod=NULL, par.prior = TRUE, prior.plots = FALSE)
    ret[["combat"]]=datCombat


    if(length(genesUse)>1)
    {
        dat=dat[map_chr(rownames(dat),function(x) strsplit(x,"_")[[1]][1]) %in% genesUse,]
        datCombat=datCombat[map_chr(rownames(datCombat),function(x) strsplit(x,"_")[[1]][1]) %in% genesUse,]
    }

 


    ##save expression data
    num=dim(dat)[2]
    dat=data.frame(dat)
    dat["geneid"]=rownames(dat)
    dat=dat[,c((num+1),1:num)]
    write.table(dat,expressFile,sep="\t",quote=F,row.names=F) 

    dat=datCombat
    num=dim(dat)[2]
    dat=data.frame(dat)
    dat["geneid"]=rownames(dat)
    dat=dat[,c((num+1),1:num)]
    write.table(dat,sub("txt$","combat.txt",expressFile),sep="\t",quote=F,row.names=F) 

    return(ret)
}



##Makes meta file
makeMeta<-function(pca,meta_PB,metaFile,covsToUse=c("sex","age",sub("^","PC",1:5)))
{
    pc=pca$eigenvect
    colnames(pc)=sub("^","PC",1:dim(pc)[2])
    rownames(pc)=pca[[1]]
    pc=pc[rownames(meta_PB),]
    comb=cbind(meta_PB,pc)
    comb=comb[,covsToUse]
    comb=model.matrix(as.formula(paste0("~",paste0(covsToUse,collapse="+"))),comb)
    comb=comb[,2:dim(comb)[2]]
    comb=data.frame(t(comb))
    #comb=comb[,c(dim(comb)[2],1:(dim(comb)[2]-1))]
    write.table(comb,metaFile,sep="\t",row.names=F,quote=F)
    return(cbind(meta_PB,pc))
}

FullPipeline_Auto<-function()
{
    meta=qread("meta.qs")
    meta=meta[meta$MajorCellTypes=="GLU_Neurons",]
    snps=qread("snps.qs")
    genes=qread("genes.qs")
    gdsFile="temp.gds"
    FullPipeline(meta,snpsUse=snps,genesUse=genes,gdsFile=gdsFile,DSOnly=T)
}

FullPipeline<-function(meta,vcfFile="/stanley/levin_asap_storage/612-eqtl/GenotypeData/Clean_vcf/Combined/comb_new.no.chr.vcf.gz",tmpDir="temp",snpsUse=c(),genesUse=c(),covsToUse=c("sex","age","RIN","PMI",sub("^","PC",1:5)),gdsFile="",numDS=-1,retFile="",DSOnly=F)
{
    print("Get GDS")
    if(nchar(gdsFile)<5)
    {
        gdsFile=paste0(tmpDir,"/temp.gds")
    }
    if(!file.exists(gdsFile))
    {
        gds=makeGDS(vcfFile,gdsFile=gdsFile)
    }
    #pca=getPCA(gdsFile)

    print("Make pseudobulk")
    expressFile=paste0(tmpDir,"/Expression.txt")
    if(nchar(retFile)<1)
    {
        retFile=paste0(tmpDir,"/Expression.qs")
    }

    if(!file.exists(expressFile))
    {
        ret=makePB(meta,expressFile,genesUse,numDS=numDS,retFile=retFile)
        qsave(ret,retFile)
    }
    ret=qread(retFile)
    
    dat=ret[[1]]
    


    print("Get Metadata")
    metaFile=paste0(tmpDir,"/Meta.txt")
    if(!file.exists(metaFile))
    {
        print("Make PCA")
        pca=getPCA(gdsFile)
        print("Save meta")
        metaPB=makeMeta(pca,ret[[2]],metaFile,covsToUse=covsToUse)
    }

    print("Get SNP info")
    snpFile=paste0(tmpDir,"/SNPs.txt")
    if(!file.exists(snpFile))
    {
        snps=makeSNPFile(vcfFile,dat,snpFile,snpsUse)
    }

    print("Run test!")
    ret=NULL
    if(!DSOnly)
    {
        ret=Run_eQTL(expressFile,snpFile,metaFile,tmpDir,snpsUse,genesUse)
        qsave(ret,paste0(tmpDir,"/all.test.qs"))
    }
    expressFile=sub("txt$","combat.txt",expressFile)
    print("DS")
    ret2=RunDS_eQTL(expressFile,snpFile,metaFile,tmpDir)
    return(ret)

}

Run_eQTL<-function(expressFile,snpFile,metaFile,tmpDir,snpsUse=c(),genesUse=c())
{

    ret=list()
    print("Prep SNPs")
    snps = SlicedData$new()
    snps$fileDelimiter = "\t"
    snps$fileOmitCharacters = "NA"
    snps$fileSkipRows = 1
    snps$fileSkipColumns = 1
    snps$fileSliceSize = 2000
    snps$LoadFile(snpFile)

    print("Prep covariance")
    gene = SlicedData$new()
    gene$fileDelimiter = "\t"
    gene$fileOmitCharacters = "NA"
    gene$fileSkipRows = 1
    gene$fileSkipColumns = 1
    gene$fileSliceSize = 2000
    gene$LoadFile(expressFile)

    print("Prep expression")
    cvrt = SlicedData$new()
    cvrt$fileDelimiter = "\t"
    cvrt$fileOmitCharacters = "NA"
    cvrt$fileSkipRows = 1
    cvrt$fileSkipColumns = 0
    cvrt$fileSliceSize = 2000
    cvrt$LoadFile(metaFile)

    useModel = modelLINEAR
    pvOutputThreshold = 1;
    errorCovariance = numeric();
    useModel = modelLINEAR
    output_file_name="temp.txt"

    print("Run!")
    me = Matrix_eQTL_engine(
        snps = snps,
        gene = gene,
        cvrt = cvrt,
        output_file_name = output_file_name,
        pvOutputThreshold = pvOutputThreshold,
        useModel = useModel,
        errorCovariance = errorCovariance,
        verbose = TRUE,
        pvalue.hist = TRUE,
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = FALSE);
    qsave(me,paste0(tmpDir,"/eQTL.CPM.Linear.qs"))
    mrk=me$all$eqtls
    mrk["Gene"]=map_chr(mrk$gene,function(x) strsplit(x,"_")[[1]][1])
    mrk["Match"]="No"
    for(i in 1:length(snpsUse))
    {
        mrk[mrk$snps==snpsUse[i] & mrk$Gene==genesUse[i],"Match"]="Yes"
    }
    mrk=mrk[mrk$Match=="Yes",]
    mrk["padj"]=p.adjust(mrk$pvalue,"fdr")
    ret[["CPM_Linear"]]=mrk

    print("Run!")
    useModel = modelANOVA
    me = Matrix_eQTL_engine(
        snps = snps,
        gene = gene,
        cvrt = cvrt,
        output_file_name = output_file_name,
        pvOutputThreshold = pvOutputThreshold,
        useModel = useModel,
        errorCovariance = errorCovariance,
        verbose = TRUE,
        pvalue.hist = TRUE,
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = FALSE);
    qsave(me,paste0(tmpDir,"/eQTL.CPM.Anova.qs"))
    mrk=me$all$eqtls
    mrk["Gene"]=map_chr(mrk$gene,function(x) strsplit(x,"_")[[1]][1])
    mrk["Match"]="No"
    for(i in 1:length(snpsUse))
    {
        mrk[mrk$snps==snpsUse[i] & mrk$Gene==genesUse[i],"Match"]="Yes"
    }
    mrk=mrk[mrk$Match=="Yes",]
    mrk["padj"]=p.adjust(mrk$pvalue,"fdr")
    ret[["CPM_Anova"]]=mrk

    print("ComBat")
    expressFile=sub("txt$","combat.txt",expressFile)
    print("Prep covariance")
    gene = SlicedData$new()
    gene$fileDelimiter = "\t"
    gene$fileOmitCharacters = "NA"
    gene$fileSkipRows = 1
    gene$fileSkipColumns = 1
    gene$fileSliceSize = 2000
    gene$LoadFile(expressFile)

    pvOutputThreshold = 1;
    errorCovariance = numeric();
    useModel = modelLINEAR
    output_file_name="temp.txt"

    print("Run!")
    me = Matrix_eQTL_engine(
        snps = snps,
        gene = gene,
        cvrt = cvrt,
        output_file_name = output_file_name,
        pvOutputThreshold = pvOutputThreshold,
        useModel = useModel,
        errorCovariance = errorCovariance,
        verbose = TRUE,
        pvalue.hist = TRUE,
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = FALSE);
    qsave(me,paste0(tmpDir,"/eQTL.ComBat.Linear.qs"))
    mrk=me$all$eqtls
    mrk["Gene"]=map_chr(mrk$gene,function(x) strsplit(x,"_")[[1]][1])
    mrk["Match"]="No"
    for(i in 1:length(snpsUse))
    {
        mrk[mrk$snps==snpsUse[i] & mrk$Gene==genesUse[i],"Match"]="Yes"
    }
    mrk=mrk[mrk$Match=="Yes",]
    mrk["padj"]=p.adjust(mrk$pvalue,"fdr")
    ret[["Combat_Linear"]]=mrk

    print("Run!")
    useModel = modelANOVA
    me = Matrix_eQTL_engine(
        snps = snps,
        gene = gene,
        cvrt = cvrt,
        output_file_name = output_file_name,
        pvOutputThreshold = pvOutputThreshold,
        useModel = useModel,
        errorCovariance = errorCovariance,
        verbose = TRUE,
        pvalue.hist = TRUE,
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = FALSE);
    qsave(me,paste0(tmpDir,"/eQTL.ComBat.Anova.qs"))
    mrk=me$all$eqtls
    mrk["Gene"]=map_chr(mrk$gene,function(x) strsplit(x,"_")[[1]][1])
    mrk["Match"]="No"
    for(i in 1:length(snpsUse))
    {
        mrk[mrk$snps==snpsUse[i] & mrk$Gene==genesUse[i],"Match"]="Yes"
    }
    mrk=mrk[mrk$Match=="Yes",]
    mrk["padj"]=p.adjust(mrk$pvalue,"fdr")
    ret[["Combat_ANOVA"]]=mrk

    return(ret)



}



RunDS_eQTL<-function(expressFile,snpFile,metaFile,tmpDir)
{
    listDS=list()
    out=map(10*(2:9),function(num)
    {
        print(num)
        dat=read.table(expressFile,header=T)
        meta=read.table(metaFile,header=T)
        snpsDat=read.table(snpFile,header=T)

        samps=sample(colnames(meta),num)
        listDS[[num]]=samps
        rownames(dat)=dat[,1]
        rownames(snpsDat)=snpsDat[,1]
        dat=dat[,samps]
        snpsDat=snpsDat[,samps]
        meta=meta[,samps]
        dat=as.matrix(dat)
        snpsDat=as.matrix(snpsDat)
        meta=as.matrix(meta)


        snps=SlicedData$new()
        snps$CreateFromMatrix( snpsDat );
        cvrt=SlicedData$new()
        cvrt$CreateFromMatrix( meta );
        gene=SlicedData$new()
        gene$CreateFromMatrix( dat );



        #useModel = modelLINEAR
        pvOutputThreshold = 1;
        errorCovariance = numeric();
        useModel = modelLINEAR
        output_file_name="temp.txt"
        useModel = modelANOVA
        me = Matrix_eQTL_engine(
            snps = snps,
            gene = gene,
            cvrt = cvrt,
            output_file_name = output_file_name,
            pvOutputThreshold = pvOutputThreshold,
            useModel = useModel,
            errorCovariance = errorCovariance,
            verbose = TRUE,
            pvalue.hist = TRUE,
            min.pv.by.genesnp = FALSE,
            noFDRsaveMemory = FALSE);
        #qsave(me,paste0(tmpDir,"/eQTL.CPM.Anova.qs"))
        ret=list()
        ret[["me"]]=me
        ret[["samp"]]=samps
        return(ret)

    })
    names(out)=10*(2:9)
    qsave(out,paste0(tmpDir,"/DS.qs"))
    #qsave(listDS,paste0(tmpDir,"/DSSamps.qs"))
    return(out)
}


LoadExpressionPB<-function(meta,cbc_col="CBC",raw_col="Raw",samp_col="batchSamp",genes_include=NULL)
{   
    print("Get samps")
    samps=unique(meta[,raw_col])
    print("Rename")
    #rownames(meta)=apply(meta,1,function(x){paste(x[samp_col],x[cbc_col],sep="_")})
    print("Next")
    out=lapply(samps,function(x){
        print(x);
        dat=LoadCountMatrix(x)
        meta_tmp=meta[meta[,raw_col]==x,]
        cbc=meta_tmp[,cbc_col]
        #nams=rownames(meta_tmp)
        #names(nams)=cbc
        print("Number in meta")
        print(dim(meta_tmp))
        print("Number in dat before:")
        print(dim(dat))
        print(head(meta_tmp))
        print(head(colnames(dat)))
        dat=dat[,cbc]
        print("Number in dat after:")
        print(dim(dat))

        #colnames(dat)=as.character(nams[colnames(dat)])
        if(!is.null(genes_include))
        {
            genes=intersect(rownames(dat),genes_include)
            dat=dat[genes,]
        }
        tab=data.frame(rowSums(dat))
        colnames(tab)[1]=x 
        return(tab)
    })


    dat=do.call(cbind,out)
    colnames(dat)=map_chr(samps,function(x){
        meta_tmp=meta[meta[,raw_col]==x,]
        meta_tmp[1,"VCF_Name"]
    })
    #dat=dat[,rownames(meta)]
    return(dat)
}



CleanUp<-function(me,snpsUse,genesUse)
{
    mrk=me$all$eqtls
    print("Gene name")
    mrk["Gene"]=map_chr(mrk$gene,function(x) strsplit(x,"_")[[1]][1])
    mrk["Match"]="No"
    print("SNP")
    for(i in 1:length(snpsUse))
    {
        mrk[mrk$snps==snpsUse[i] & mrk$Gene==genesUse[i],"Match"]="Yes"
    }
    mrk=mrk[mrk$Match=="Yes",]
    mrk["FDR"]=p.adjust(mrk$pvalue,"fdr")
    return(mrk)

}