#' @title A R-based script to compute the statical significance using MR.
#' @author Hongwei Wang <\email{whweve@163.com}> =================================== load library ============================================= 
library(data.table)
library(sqldf)
library(speedglm)
library(RcppArmadillo)


# =================================== load library ============================================= 
RSQUARE <- function(object) {
  f <- object$fitted.values    
  r <- object$residuals         
  mss <- if (object$intercept) sum((f - mean(f))^2) else sum(f^2)     
  rss <- sum(r^2)   
  r.squared <- mss/(mss + rss) 
  
  #df.int <- if (object$intercept) 1L else 0L    
  #n <- length(f)     
  #rdf <- object$df     
  #adj.r.squared <- 1 - (1 - r.squared) * ((n - df.int)/rdf)
  return(r.squared)
}
# =================================== meta table   ============================================= 
#read in the metadata

metatable <- read.table(text="
                        treat localtion type name
                        ww D:/whw/qfdata/FPKM_data FPKM ww_genes_FPKM_197.txt
                        70 D:/whw/qfdata/FPKM_data FPKM 70_FPKM_186.txt
                        58 D:/whw/qfdata/FPKM_data FPKM 58_FPKM_189.txt
                        ww D:/whw/qfdata/pheno pheno WW_50_0.05_normalizeQuantiles.qqnorm
                        70 D:/whw/qfdata/pheno pheno 70_50_0.05_normalizeQuantiles.qqnorm
                        58 D:/whw/qfdata/pheno pheno 58_50_0.05_normalizeQuantiles.qqnorm
                        ww D:/whw/qfdata/topeQTL topeqtl ww.csv
                        70 D:/whw/qfdata/topeQTL topeqtl 70.csv
                        58 D:/whw/qfdata/topeQTL topeqtl 58.csv
                        ww D:/whw/qfdata/genotypedata snp WW_eGWAS.genotype
                        70 D:/whw/qfdata/genotypedata snp 70-70R_eGWAS.genotype
                        58 D:/whw/qfdata/genotypedata snp 58-58R_eGWAS.genotype
                        ",header=T)
metatable$file <- paste0(metatable$localtion,"/",metatable$name)

# =================================== define treat   ============================================ 
for (treat in c("ww","70","58")) {
  # =================================== read data      ============================================ 

eqtlfile <- metatable[metatable$treat == treat & metatable$type=="topeqtl","file"]
expfile <- metatable[metatable$treat == treat & metatable$type=="pheno","file"]
snpfile <- metatable[metatable$treat == treat & metatable$type=="snp","file"]
fpkmfile <- metatable[metatable$treat == treat & metatable$type=="FPKM","file"]
#read in the data
eqtl <- fread(eqtlfile,header=T)
#filter the top-eqtl
eqtl <- as.data.frame(eqtl)
topeqtl <- sqldf("select Gene,eQTL,SNP,min(pvalue) minpvalue from eqtl  group by Gene")

#read in the expression data
exp <- fread(expfile,header=T)
exp <- as.data.frame(exp)

pheno <- read.table("D:\\whw\\fromD\\gemma\\gwa\\file\\pheno2-r.txt",header=T)
pheno <- pheno[,c("inbred","BLUP")]
#read in the snp data
snp <- fread(snpfile,header=T)
snp <- as.data.frame(snp)

#check the snp whether is bi-allele
range(nchar(snp$allele)) #it turned out to be 3:3, all of them are bi-allele

#convert the SNP to numeric format
major_allele <- paste0(substr(snp$allele,1,1),substr(snp$allele,1,1))
#if allele equal to major allele, 0, else 2
for (j in 12:dim(snp)[2]) {
  major_allele_pos <- (snp[,j] == major_allele)
  snp[major_allele_pos,j] = 2
  snp[!major_allele_pos,j] = 0
  snp[,j] <- as.numeric(snp[,j])
}

snp <- snp[,-c(2:11)]
# =================================== calculate      =========================================== 

hold <- data.frame()

for (i in 1:dim(topeqtl)[1]) {
  #get the gene names
  gene <- topeqtl[i,1]
  #get the gene-leadsnp infomation
  gene_topeqtl_snp <- topeqtl[i,3]
  #get the exp data
  gene_pos_in_exp <- which(names(exp) == gene)
  
  #converthe SNP, n*1 format
  tmp_snp <- t(snp[snp$rs == gene_topeqtl_snp,])
  #to go on to the next if the gene has no expression or SNP data
  if (length(gene_pos_in_exp) == 0 ) {
    hold_tmp <- data.frame(gene=gene,
                                  gene_topeqtl_snp=gene_topeqtl_snp,
                                  phenotype_expression="NO-EXPRESSION",
                                  TMR="NONE",
                                  pvalue="NONE")
    hold <- rbind(hold,hold_tmp)
    next
  } else if (length(which(snp$rs == gene_topeqtl_snp)) == 0) {
    hold_tmp <- data.frame(gene=gene,
                           gene_topeqtl_snp=gene_topeqtl_snp,
                           TMR="NONE",
                           phenotype_expression="NO-SNP",
                           pvalue="NONE")
    hold <- rbind(hold,hold_tmp)
    next
  }   
    tmp_snp <- data.frame(inbred = row.names(tmp_snp),
                          #convert to numeric to avoid unexpected transform
                          genotype = as.numeric(as.character(tmp_snp[,1])))
    
    #convert the exp, n*1 format
    tmp_exp <- data.frame(inbred = exp[,1],
                          expression = exp[,gene_pos_in_exp])
    
    #construct the design matrix
    design_matrix <- merge(tmp_snp,tmp_exp,by = "inbred")
    design_matrix <- merge(design_matrix,pheno,by="inbred")
    design_matrix$intercept <- rep(1,dim(design_matrix)[1])

        #compute the effect
    #estimates of expression on genotype
    
    #lm_genotype_expression = ols.fit1(y=design_matrix$expression,x=cbind(design_matrix$intercept,design_matrix$genotype))
    lm_genotype_expression = lm(expression~genotype,data=design_matrix)
    effect_genotype_expression = lm_genotype_expression$coefficients[2]
    #estimate of phenotype on genotype
    #lm_genotype_phenotype = ols.fit1(y=design_matrix$BLUP,x=cbind(design_matrix$intercept,design_matrix$genotype))
    lm_genotype_phenotype = lm(BLUP ~ genotype,data = design_matrix)
    effect_genotype_phenotype = lm_genotype_phenotype$coefficients[2]
    #efficiency_phenotype_on_expression
    phenotype_expression = effect_genotype_phenotype/effect_genotype_expression
    
    #estimate the variance
    #lm_expression_phenotype <- lm.fit(y=design_matrix$BLUP,x=cbind(design_matrix$intercept,design_matrix$expression))
    lm_expression_phenotype <- lm(BLUP ~ expression,data = design_matrix)
    rs_expression_phenotype = summary(lm_expression_phenotype)$r.square
    rs_genotype_expression = summary(lm_genotype_expression)$r.square
    #rs_expression_phenotype = RSQUARE(lm_expression_phenotype)
    #rs_genotype_expression = RSQUARE(lm_genotype_expression)
    
  
    upper = var(design_matrix$BLUP)*(1-rs_expression_phenotype)
    down = dim(design_matrix)[1]*var(design_matrix$expression)*rs_genotype_expression
    var_phenotype_expression = upper/down
    
    TMR = phenotype_expression^2/var_phenotype_expression
    
    pvalue = pchisq(TMR,1,lower.tail = F)
    
    hold_tmp <- data.frame(gene=gene,
                           gene_topeqtl_snp=gene_topeqtl_snp,
                           phenotype_expression=phenotype_expression,
                           TMR=TMR,
                           pvalue=pvalue)
    hold <- rbind(hold,hold_tmp)
}

write.csv(hold,paste0("D:\\whw\\Rpro\\bj_wgcna\\",treat,".csv"),quote=F,row.names=F)

}


