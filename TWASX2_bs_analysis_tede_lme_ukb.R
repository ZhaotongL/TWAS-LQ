library(games)
library(nlme)
library(sandwich)
library(TEDE)
library(aSPU)
library(plink2R)
library(dplyr)
library(tidyr)
library(data.table)
source("allele_qc.R")
source("my_TEDE.R")

chr = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(chr)
chr0 = paste0('chr',chr)
### ukb phenotype
ukb_pheno_all = fread('/home/panwei/lin00374/TWAS_XSquare/ukb/ADNI_Stage1/ukb_hdl.txt')
ukb_pheno_all = na.omit(ukb_pheno_all)


###UKB
bed_dir = '/home/panwei/shared/UKBiobankIndiv/genetic_data/ukb_cal_chr'
prefix = paste0(bed_dir,chr,'_v2')

### Covariates
cov5table = read.table("/home/panwei/lin00374/TWAS_XSquare/ukb/gtex_v8_fusion/Whole_Blood.v8.covariates.txt",header = T)
### Gene Expression
gene_exp = fread("/home/panwei/shared/GTEx-v8/GTEx_Analysis_v8_expression/expression_matrices/Whole_Blood.v8.normalized_expression.bed.gz")
colnames(gene_exp)[1] = 'chr'
gene_exp = gene_exp %>% filter(chr==chr0)
write.table(gene_exp$gene_id,paste0('./result/gene_',chr0,'.txt'),row.names=F,quote=F)
#write.table(cbind(0,colnames(gene_exp)[-(1:4)]),'/home/panwei/lin00374/TWAS_XSquare/ukb/gtex_v8_fusion/Whole_Blood_buddy.txt',quote=F,row.names=F,col.names=F)
whole_blood_buddy = '/home/panwei/lin00374/TWAS_XSquare/ukb/gtex_v8_fusion/Whole_Blood_buddy.txt'

cat(nrow(gene_exp),"\n")
real_data_result = list()
# start analysis ----------------------------------------------------------
for(gene_ind in 1:nrow(gene_exp)){
    
  cat(gene_ind,"\n")
  ukb_pheno = ukb_pheno_all
  gene_name = as.character(gene_exp$gene_id[gene_ind])
  start = floor(gene_exp$start[gene_ind]/1000) - 100
  end = ceiling(gene_exp$end[gene_ind]/1000) + 100
  
  plink_command = paste0("module load plink; \n","plink --bfile /home/panwei/lin00374/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased_MAF01",
                        " --chr ",chr," --from-kb ",
                        start," --to-kb ",end,
                        " --keep ",whole_blood_buddy,
                        " --geno 0 --maf 0.05 --hwe 0.001  ",
                        " --make-bed --out For_Individual_Genes/",
                        gene_name,sep="")
  plink_msg=system(plink_command,intern=TRUE)
  
  snp = tryCatch(read_plink(paste("For_Individual_Genes/",gene_name,sep="")),
                 error=function(e){cat("ERROR :",
                                       conditionMessage(e),
                                       "\n")})

  remove_command = paste("rm For_Individual_Genes/",
                         gene_name,".*",sep="")
  system(remove_command)

  if(is.null(snp))
  {
    FINAL_RESULT = NULL
    real_data_result = c(real_data_result,list(FINAL_RESULT = FINAL_RESULT))
    
    next()
  }
  write.table(snp$bim$V2,paste0('ForUKB/',gene_name,'_rs.txt'),row.names=F,quote=F,col.names=F)

  grep_command = paste0("grep -F -m ",length(snp$bim$V2)," -f ForUKB/",gene_name,"_rs.txt /home/panwei/shared/GTEx-v8/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt  > ForUKB/",gene_name,"_rs1.txt")
  system(grep_command)
  gtex_dic = read.table(paste0('ForUKB/',gene_name,'_rs1.txt'))
  write.table(gtex_dic$V7,paste0('ForUKB/',gene_name,'_rs.txt'),row.names=F,quote=F,col.names=F)

  ukb_out = paste0('ForUKB/',gene_name)
  ukb_plink_command = paste0('module load plink; plink --bfile ',prefix,' --chr ',chr,
                          ' --keep /home/panwei/lin00374/TWAS_XSquare/ukb/white_unrelated_keep_ind.txt ',
                          ' --extract ','ForUKB/',gene_name,'_rs.txt',
                          ' --mind 0 --make-bed',
                          ' --out ',ukb_out)
  ukb_plink_msg = system(ukb_plink_command,intern=TRUE)

  ukb_snp = tryCatch(read_plink(paste("ForUKB/",gene_name,sep="")),
                 error=function(e){cat("ERROR :",
                                       conditionMessage(e),
                                       "\n")})
  
  remove_command = paste("rm ForUKB/",
                         gene_name,".*",sep="")
  system(remove_command)
  remove_command = paste("rm ForUKB/",
                         gene_name,"_rs*",sep="")
  system(remove_command)

  if(is.null(ukb_snp))
  {
    FINAL_RESULT = NULL
    real_data_result = c(real_data_result,list(FINAL_RESULT = FINAL_RESULT))

    next()
  }



  
  ### overlap between ADNI SNP and UKB
#  colnames(snp$bed) = snp$bim$V4
  rownames(ukb_snp$bed) = ukb_snp$fam[,2]
  
#  colnames(ukb_snp$bim)[4:6] = c("Position","A1","A2")
  gtex_bim_tmp = merge(snp$bim,gtex_dic,by.x='V2',by.y='V1',sort=F)
  snp$bim[,2] = gtex_bim_tmp$V7
  colnames(snp$bed) = snp$bim[,2]
  rownames(snp$bed) = snp$fam[,2]
  snp_bim_ukb = merge(snp$bim,ukb_snp$bim,by=c("V1","V2"))
  remove_flip = allele.qc(snp_bim_ukb$V5.x,snp_bim_ukb$V6.x,
                          snp_bim_ukb$V5.y,snp_bim_ukb$V6.y)
  flip_snp = snp_bim_ukb$V2[remove_flip$flip]
  snp_bim_ukb = snp_bim_ukb[remove_flip$keep,]

  ukb_snp_bed = ukb_snp$bed[,which(colnames(ukb_snp$bed)%in% snp_bim_ukb$V2),drop=FALSE] 
  ukb_snp_bed[,which(colnames(ukb_snp_bed) %in% flip_snp)] = 2 - ukb_snp_bed[,which(colnames(ukb_snp_bed) %in% flip_snp)]
  
  snp_bed_ind = NULL
  snp_bed_colname = colnames(snp$bed)
  for(i in 1:nrow(snp_bim_ukb))
  {
    snp_bed_ind = c(snp_bed_ind,which(snp_bed_colname == colnames(ukb_snp_bed)[i]))
  }
  SNP_BED = snp$bed[,snp_bed_ind,drop=FALSE]
  cov5table = cov5table[match(rownames(SNP_BED),cov5table$ID),]


  keep_ukb_indiv = intersect(ukb_pheno$f.eid,ukb_snp$fam[,2])
  ukb_pheno = ukb_pheno %>% filter(f.eid %in% keep_ukb_indiv)
  ukb_snp_bed = ukb_snp_bed[which(rownames(ukb_snp_bed)%in%keep_ukb_indiv),,drop=FALSE]
  ukb_pheno = ukb_pheno[match(rownames(ukb_snp_bed),ukb_pheno$f.eid),]

  if(nrow(ukb_pheno)< 3){
      FINAL_RESULT = NULL
      real_data_result = c(real_data_result,list(FINAL_RESULT = FINAL_RESULT))
      next()
  }
  
  
  if(length(snp_bed_ind) < 2)
  {
    FINAL_RESULT = NULL
    real_data_result = c(real_data_result,list(FINAL_RESULT = FINAL_RESULT))
    
    next()
  }
  ### prune
  cor_cutoff = 0.8
  cor_bed = abs(cor(SNP_BED))
  cor_bed = (cor_bed < cor_cutoff)^2
  diag(cor_bed) = 1
  i = 1
  while(i < nrow(cor_bed) )
  {
    ind = which(cor_bed[i,] == 1)
    cor_bed = as.matrix(cor_bed[ind,ind])
    i = i + 1
  }
  if(nrow(cor_bed) == 1)
  {
    FINAL_RESULT = NULL
    real_data_result = c(real_data_result,list(FINAL_RESULT = FINAL_RESULT))
    
    next()
  }
  ind = which(is.element(colnames(SNP_BED),colnames(cor_bed)))
  SNP_BED = SNP_BED[,ind,drop=F]
  ukb_snp_bed = ukb_snp_bed[,ind]
  snp_bim_ukb = snp_bim_ukb[ind,,drop=F]

  
  ### regress gene exp on covariates
  
  X = as.matrix(gene_exp[gene_ind,-(1:4)])
  X = X[,match(rownames(SNP_BED),colnames(X)),drop=F]
  X = scale(as.numeric(X))
#  X = as.numeric(gene_exp[gene_ind,-(1:4)]) #A1,A3
  X1 = X
  X2 = X^2 #A0
  
  lm_RegCov = lm(X1 ~.,data = data.frame(X1,cov5table[,-c(1,2)]))
  X1 = lm_RegCov$residuals
  
#A0
  lm_RegCov = lm(X2 ~.,data = data.frame(X2,cov5table[,-c(1,2)]))
  X2 = lm_RegCov$residuals
##
    
  X1 = scale(X1) #A0,A1
  X2 = scale(X2) #A0
#  X2 = X1^2 #A1,A2,A3
  
  if(ncol(SNP_BED) > 50)
  {
    ind1 = order(abs(cor(X1,SNP_BED)),decreasing = T)[1:50]
    ind2 = order(abs(cor(X2,SNP_BED)),decreasing = T)[1:50]
    ind = union(ind1,ind2)
    
    ind = sort(ind)
    SNP_BED = SNP_BED[,ind]
    ukb_snp_bed = ukb_snp_bed[,ind]
    snp_bim_ukb = snp_bim_ukb[ind,]
  }
  
  SNP_BED = scale(SNP_BED)
  ukb_snp_bed = scale(ukb_snp_bed)

  
  ###
  
  lm_stage1_X1 = lm(X1 ~ SNP_BED)
 # lm_stage1_X2 = lm(X2 ~ SNP_BED)
  
  hatbetaX1 = lm_stage1_X1$coefficients[-1]
  
  na_ind = which(!is.na(hatbetaX1))
  check_ukb = apply(ukb_snp_bed,2,sd)
  na_ind_ukb = which(!is.na(check_ukb))
  na_ind = intersect(na_ind,na_ind_ukb)
  SNP_BED = SNP_BED[,na_ind,drop=FALSE]
  ukb_snp_bed = ukb_snp_bed[,na_ind,drop=FALSE]

  
  lm_stage1_X1 = 
    step(lm(X1~.,data = data.frame(X1,SNP_BED)),direction = "backward",trace=FALSE)
  AIC_stage1_X1 = AIC(lm_stage1_X1)
  BIC_stage1_X1 = BIC(lm_stage1_X1)
  lm_stage1_X1 = summary(lm_stage1_X1)
  rsq_stage1_X1 = lm_stage1_X1$r.squared
  adjrsq_stage1_X1 = lm_stage1_X1$adj.r.squared
  se_betaX1 = hatbetaX1 = rep(0,ncol(SNP_BED))
  coef_X1 = lm_stage1_X1$coefficients
  name_coef_X1 = substr(rownames(coef_X1),1,30)
  name_SNP_BED = colnames(SNP_BED)
  for(beta_ind in 1:nrow(coef_X1))
  {
    ii = which(name_SNP_BED == name_coef_X1[beta_ind])
    hatbetaX1[ii] = coef_X1[beta_ind,1]
    se_betaX1[ii] = coef_X1[beta_ind,2]
  }
  

  lm_stage1_X2 = 
    step(lm(X2~.,data = data.frame(X2,SNP_BED)),direction = "backward",trace=FALSE)
  AIC_stage1_X2 = AIC(lm_stage1_X2)
  BIC_stage1_X2 = BIC(lm_stage1_X2)
  lm_stage1_X2 = summary(lm_stage1_X2)
  rsq_stage1_X2 = lm_stage1_X2$r.squared
  adjrsq_stage1_X2 = lm_stage1_X2$adj.r.squared
  se_betaX2 = hatbetaX2 = rep(0,ncol(SNP_BED))
  coef_X2 = lm_stage1_X2$coefficients
  name_coef_X2 = substr(rownames(coef_X2),1,30)
  name_SNP_BED = colnames(SNP_BED)
  for(beta_ind in 1:nrow(coef_X2))
  {
    ii = which(name_SNP_BED == name_coef_X2[beta_ind])
    hatbetaX2[ii] = coef_X2[beta_ind,1]
    se_betaX2[ii] = coef_X2[beta_ind,2]
  }
  

  if((sum(abs(hatbetaX1)) + sum(abs(hatbetaX2))) == 0 )
  {
    FINAL_RESULT = NULL
    real_data_result = c(real_data_result,list(FINAL_RESULT = FINAL_RESULT))
    
    next()
  }

### UKB 
  Xhat_ukb = ukb_snp_bed %*% hatbetaX1
  Xhat2_ukb = Xhat_ukb^2
  X2hat_ukb = ukb_snp_bed %*% hatbetaX2

  y_ukb = scale(ukb_pheno$hdl)
  n_ukb = length(y_ukb)
  YY = crossprod(y_ukb)
  Id <- factor(rep(1, length = n_ukb))

  
##M1
  X1snps = which(hatbetaX1!=0)
  stage2X1_lmm = lm(y_ukb~Xhat_ukb)
  stage2X1_lm = summary(stage2X1_lmm)
  if(length(X1snps)>1){
  GY1 = crossprod(ukb_snp_bed[,X1snps,drop=FALSE],y_ukb)
  LDcov1 = crossprod(ukb_snp_bed[,X1snps,drop=FALSE])
  mu1 = stage2X1_lmm$fitted.values
  Gmu1 = t(ukb_snp_bed[,X1snps,drop=FALSE]) %*% mu1
  effect_GX1 = hatbetaX1[X1snps]
  my_TEDE_Sc1 = my_TEDE_Sc(effect_GX=effect_GX1,
                    Gmu=Gmu1,G=ukb_snp_bed[,X1snps,drop=FALSE],
                    LDcov=LDcov1,GY=GY1,YY=YY,N=n_ukb,rcx1=vcov(lm_stage1_X1)[-1,-1])
  my_TEDE_aSPU1 = my_TEDE_aSPU(effect_GX=effect_GX1,
                    Gmu=Gmu1,G=ukb_snp_bed[,X1snps,drop=FALSE],
                    LDcov=LDcov1,GY=GY1,YY=YY,N=n_ukb,rcx1=vcov(lm_stage1_X1)[-1,-1])
  X1snps_mat = ukb_snp_bed[,X1snps,drop=FALSE]
  X1snps_sum_vec = apply(X1snps_mat,1,sum)
  m1.block<-list(Id = pdIdent(~X1snps_mat-1))
  stage2X1_egger_lme = tryCatch(lme(y_ukb ~ Xhat_ukb + X1snps_sum_vec, random = m1.block,method='ML',control = lmeControl(opt = "optim")),
                    error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
  stage2X1_ivw_lme = tryCatch(lme(y_ukb ~ Xhat_ukb , random = m1.block,method='ML',control = lmeControl(opt = "optim")),
                    error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
  if(is.null(stage2X1_egger_lme)){
      M1_egger_lme=M1_egger_lm_pval=NULL
  }else{
  M1_egger_lme = summary(stage2X1_egger_lme)
  stage2X1_egger_lmm = lm(y_ukb ~ Xhat_ukb + X1snps_sum_vec)
  M1_HC = vcovHC(stage2X1_egger_lmm,type='HC0')
  M1_Xhat_pval = pchisq(stage2X1_egger_lmm$coefficients[2]^2/M1_HC[2,2],df=1,lower.tail=F)
  M1_Z_pval = pchisq(stage2X1_egger_lmm$coefficients[3]^2/M1_HC[3,3],df=1,lower.tail=F)
  M1_egger_lm_pval = c(M1_Xhat_pval,M1_Z_pval)
  }
  if(is.null(stage2X1_ivw_lme)){
      M1_ivw_lme=NULL
  }else{
      M1_ivw_lme = summary(stage2X1_ivw_lme)
  }
  }else{
      my_TEDE_Sc1=my_TEDE_aSPU1=M1_egger_lme = M1_egger_lm_pval=M1_ivw_lme=NULL
  }

  stage2X1 = list(theta=stage2X1_lm$coefficients[-1,1],
                  theta_var = stage2X1_lm$coefficients[-1,2]^2,
                  theta_pval = stage2X1_lm$coefficients[-1,4],
                  sigma = stage2X1_lm$sigma,
                  fstat = stage2X1_lm$fstatistic,
                  aic = AIC(stage2X1_lmm),
                  bic = BIC(stage2X1_lmm),
                  my_TEDE_Sc=my_TEDE_Sc1,
                  my_TEDE_aSPU=my_TEDE_aSPU1,
                  ivw_lme_pval = M1_ivw_lme$tTable[2,5],
                  egger_lm_pval = M1_egger_lm_pval,
                  egger_lme_pval = M1_egger_lme$tTable[-1,5]
                   ) 
                 
#M2.1
  X2snps = which(hatbetaX2!=0)
  stage2X2_lmm = lm(y_ukb~X2hat_ukb) 
  stage2X2_lm = summary(stage2X2_lmm)
  if(length(X2snps)>1){
  GY2 = crossprod(ukb_snp_bed[,X2snps,drop=FALSE],y_ukb)
  LDcov2 = crossprod(ukb_snp_bed[,X2snps,drop=FALSE])
  effect_GX2 = hatbetaX2[X2snps]
  mu2 = stage2X2_lmm$fitted.values
  Gmu2 = crossprod(ukb_snp_bed[,X2snps,drop=FALSE],mu2)
  my_TEDE_Sc2.1 = my_TEDE_Sc(effect_GX=effect_GX2,
                    Gmu=Gmu2,
                    LDcov=LDcov2,GY=GY2,YY=YY,N=n_ukb,rcx1=vcov(lm_stage1_X2)[-1,-1])
  my_TEDE_aSPU2.1 = my_TEDE_aSPU(effect_GX=effect_GX2,
                    Gmu=Gmu2,
                    LDcov=LDcov2,GY=GY2,YY=YY,N=n_ukb,rcx1=vcov(lm_stage1_X2)[-1,-1])
  X2snps_mat = ukb_snp_bed[,X2snps,drop=FALSE]
      X2snps_sum_vec = apply(X2snps_mat,1,sum)
      m2.1.block<-list(Id = pdIdent(~X2snps_mat-1))
      stage2X2_egger_lme = tryCatch(lme(y_ukb ~ X2hat_ukb + X2snps_sum_vec, random = m2.1.block,control = lmeControl(opt = "optim"),method='ML'),
                            error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
      stage2X2_ivw_lme = tryCatch(lme(y_ukb ~ X2hat_ukb , random = m2.1.block,control = lmeControl(opt = "optim"),method='ML'),
                            error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
      if(is.null(stage2X2_egger_lme)){
          M2.1_egger_lme=M2.1_egger_lm_pval=NULL
      }else{
      M2.1_egger_lme = summary(stage2X2_egger_lme)
      stage2X2_egger_lmm = lm(y_ukb ~ X2hat_ukb + X2snps_sum_vec)
      M2.1_HC = vcovHC(stage2X2_egger_lmm,type='HC0')
      M2.1_Xhat_pval = pchisq(stage2X2_egger_lmm$coefficients[2]^2/M2.1_HC[2,2],df=1,lower.tail=F)
      M2.1_Z_pval = pchisq(stage2X2_egger_lmm$coefficients[3]^2/M2.1_HC[3,3],df=1,lower.tail=F)
      M2.1_egger_lm_pval = c(M2.1_Xhat_pval,M2.1_Z_pval)
    }
    if(is.null(stage2X2_ivw_lme)){
        M2.1_ivw_lme=NULL
    }else{
        M2.1_ivw_lme=summary(stage2X2_ivw_lme)
    }
  }else{
      my_TEDE_Sc2.1=my_TEDE_aSPU2.1=M2.1_egger_lme = M2.1_egger_lm_pval= M2.1_ivw_lme=NULL
  }
  stage2X2 = list(theta=stage2X2_lm$coefficients[-1,1],
                  theta_var = stage2X2_lm$coefficients[-1,2]^2,
                  theta_pval = stage2X2_lm$coefficients[-1,4],
                  sigma = stage2X2_lm$sigma,
                  fstat = stage2X2_lm$fstatistic,
                  aic = AIC(stage2X2_lmm),
                  bic = BIC(stage2X2_lmm), 
                  my_TEDE_Sc=my_TEDE_Sc2.1,
                  my_TEDE_aSPU=my_TEDE_aSPU2.1,
                  ivw_lme_pval = M2.1_ivw_lme$tTable[2,5],
                  egger_lm_pval = M2.1_egger_lm_pval,
                  egger_lme_pval = M2.1_egger_lme$tTable[-1,5]
                  )
##M3.1 
  X3snps = union(X1snps,X2snps)
  LDcov3 = crossprod(ukb_snp_bed[,X3snps,drop=FALSE])
  effect_GX3 = matrix(0,nrow=length(X3snps),ncol=2)
  ii1 = match(name_coef_X1[-1],colnames(LDcov3))
  ii2 = match(name_coef_X2[-1],colnames(LDcov3))
  effect_GX3[ii1,1] = coef_X1[-1,1]
  effect_GX3[ii2,2] = coef_X2[-1,1]
  stage2Xboth_lmm = lm(y_ukb~Xhat_ukb+X2hat_ukb) 
  stage2Xboth_lm = summary(stage2Xboth_lmm)
  if(length(X3snps)>1){
  mu3 = stage2Xboth_lmm$fitted.values
  Gmu3 = crossprod(ukb_snp_bed[,X3snps,drop=FALSE],mu3)
  GY3 = crossprod(ukb_snp_bed[,X3snps,drop=FALSE],y_ukb)
  GG1 = crossprod(ukb_snp_bed[,X3snps,drop=FALSE],ukb_snp_bed[,X1snps,drop=FALSE])
  GG2 = crossprod(ukb_snp_bed[,X3snps,drop=FALSE],ukb_snp_bed[,X2snps,drop=FALSE])
  my_TEDE_Sc3.1 = my_TEDE_Sc(effect_GX=effect_GX3,
                    Gmu=Gmu3,
                    LDcov=LDcov3,GY=GY3,YY=YY,N=n_ukb,
                    rcx1=vcov(lm_stage1_X1)[-1,-1],rcx2=vcov(lm_stage1_X2)[-1,-1],
                    GG1=GG1,GG2=GG2)
  my_TEDE_aSPU3.1 = my_TEDE_aSPU(effect_GX=effect_GX3,
                    Gmu=Gmu3,
                    LDcov=LDcov3,GY=GY3,YY=YY,N=n_ukb,
                    rcx1=vcov(lm_stage1_X1)[-1,-1],rcx2=vcov(lm_stage1_X2)[-1,-1],
                    GG1=GG1,GG2=GG2)
      X3snps_mat = ukb_snp_bed[,X3snps]
      X3snps_sum_vec = apply(X3snps_mat,1,sum)
      m3.1.block<-list(Id = pdIdent(~X3snps_mat-1))
      stage2Xboth_egger_lme = tryCatch(lme(y_ukb ~ Xhat_ukb + X2hat_ukb + X3snps_sum_vec, random = m3.1.block,control = lmeControl(opt = "optim"),method='ML'),
                            error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
      stage2Xboth_ivw_lme = tryCatch(lme(y_ukb ~ Xhat_ukb + X2hat_ukb , random = m3.1.block,control = lmeControl(opt = "optim"),method='ML'),
                            error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
      if(is.null(stage2Xboth_egger_lme)){
          M3.1_egger_lme_pval = M3.1_egger_lm_pval = NULL
      }else{
      M3.1_egger_lme = summary(stage2Xboth_egger_lme)
      M3.1_egger_lme_chisq = pchisq(t(M3.1_egger_lme$coefficients$fixed[2:3]) %*% ginv(vcov(stage2Xboth_egger_lme)[2:3,2:3]) %*% M3.1_egger_lme$coefficients$fixed[2:3],
                            df =2 ,lower.tail=F)
      M3.1_egger_lme_pval = c(M3.1_egger_lme$tTable[-1,5],M3.1_egger_lme_chisq)

      stage2Xboth_egger_lmm = lm(y_ukb ~ Xhat_ukb + X2hat_ukb + X3snps_sum_vec)
      M3.1_HC = vcovHC(stage2Xboth_egger_lmm,type='HC0')
      M3.1_Z_pval = pchisq(stage2Xboth_egger_lmm$coefficients[4]^2/M3.1_HC[4,4],df=1,lower.tail=F)
      M3.1_Xhat_pval = pchisq(t(stage2Xboth_egger_lmm$coefficients[2:3]) %*% ginv(M3.1_HC[2:3,2:3]) %*% stage2Xboth_egger_lmm$coefficients[2:3],
                                df = 2, lower.tail=F)
      M3.1_egger_lm_pval = c(M3.1_Xhat_pval,M3.1_Z_pval)
      }
      if(is.null(stage2Xboth_ivw_lme)){
          M3.1_ivw_lme_pval=NULL
      }else{
          M3.1_ivw_lme = summary(stage2Xboth_ivw_lme)
          M3.1_ivw_lme_chisq = pchisq(t(M3.1_ivw_lme$coefficients$fixed[2:3]) %*% ginv(vcov(stage2Xboth_ivw_lme)[2:3,2:3]) %*% M3.1_ivw_lme$coefficients$fixed[2:3],
                                  df =2 ,lower.tail=F)
          M3.1_ivw_lme_pval = c(M3.1_ivw_lme$tTable[-1,5],M3.1_ivw_lme_chisq)
      }
    }else{
      my_TEDE_Sc3.1=my_TEDE_aSPU3.1=M3.1_egger_lme_pval = M3.1_egger_lm_pval= M3.1_ivw_lme_pval=NULL
      }
  stage2Xboth = list(theta=stage2Xboth_lm$coefficients[-1,1],
                  theta_var = vcov(stage2Xboth_lm)[-1,-1],
                  theta_pval = stage2Xboth_lm$coefficients[-1,4],
                  sigma = stage2Xboth_lm$sigma,
                  fstat = stage2Xboth_lm$fstatistic,
                  aic = AIC(stage2Xboth_lmm),
                  bic = BIC(stage2Xboth_lmm), 
                  my_TEDE_Sc=my_TEDE_Sc3.1,
                  my_TEDE_aSPU=my_TEDE_aSPU3.1,
                  ivw_lme_pval = M3.1_ivw_lme_pval,
                  egger_lm_pval = M3.1_egger_lm_pval,
                  egger_lme_pval = M3.1_egger_lme_pval)

##M2.2
  stage2X2_lmm1 = lm(y_ukb~Xhat2_ukb)
  stage2X2_lm1 = summary(stage2X2_lmm1)
  if(length(X1snps)>1){
  mu2.2 = stage2X2_lmm1$fitted.values
  Gmu2.2 = t(ukb_snp_bed[,X1snps,drop=FALSE]) %*% mu2.2
  my_TEDE_Sc2.2 = my_TEDE_Sc(
                    Gmu=Gmu2.2,
                    LDcov=LDcov1,GY=GY1,YY=YY,N=n_ukb,sigY2=stage2X2_lm1$sigma)
  my_TEDE_aSPU2.2 = my_TEDE_aSPU(
                    Gmu=Gmu2.2,
                    LDcov=LDcov1,GY=GY1,YY=YY,N=n_ukb,sigY2=stage2X2_lm1$sigma)
  X1snps_mat = ukb_snp_bed[,X1snps,drop=FALSE]
  X1snps_sum_vec = apply(X1snps_mat,1,sum)
  m1.block<-list(Id = pdIdent(~X1snps_mat-1))
  stage2X2_egger_lme1 = tryCatch(lme(y_ukb ~ Xhat2_ukb + X1snps_sum_vec, random = m1.block,method='ML',control = lmeControl(opt = "optim")),
                    error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
  stage2X2_ivw_lme1 = tryCatch(lme(y_ukb ~ Xhat2_ukb , random = m1.block,method='ML',control = lmeControl(opt = "optim")),
                    error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
  if(is.null(stage2X2_egger_lme1)){
      M2.2_egger_lme=M2.2_egger_lm_pval=NULL
  }else{
  M2.2_egger_lme = summary(stage2X2_egger_lme1)
  stage2X2_egger_lmm1 = lm(y_ukb ~ Xhat2_ukb + X1snps_sum_vec)
      M2.2_HC = vcovHC(stage2X2_egger_lmm1,type='HC0')
      M2.2_Xhat_pval = pchisq(stage2X2_egger_lmm1$coefficients[2]^2/M2.2_HC[2,2],df=1,lower.tail=F)
      M2.2_Z_pval = pchisq(stage2X2_egger_lmm1$coefficients[3]^2/M2.2_HC[3,3],df=1,lower.tail=F)
      M2.2_egger_lm_pval = c(M2.2_Xhat_pval,M2.2_Z_pval)
  }
  if(is.null(stage2X2_ivw_lme1)){
      M2.2_ivw_lme=NULL
  }else{
      M2.2_ivw_lme = summary(stage2X2_ivw_lme1)
  }
  }else{
     my_TEDE_Sc2.2 =  my_TEDE_aSPU2.2=M2.2_egger_lme=M2.2_egger_lm_pval=M2.2_ivw_lme=NULL
  }

  stage2X2_1 = list(theta=stage2X2_lm1$coefficients[-1,1],
                  theta_var = stage2X2_lm1$coefficients[-1,2]^2,
                  theta_pval = stage2X2_lm1$coefficients[-1,4],
                  sigma = stage2X2_lm1$sigma,
                  fstat = stage2X2_lm1$fstatistic,
                  aic = AIC(stage2X2_lmm1),
                  bic = BIC(stage2X2_lmm1),
                  my_TEDE_Sc=my_TEDE_Sc2.2,
                  my_TEDE_aSPU=my_TEDE_aSPU2.2,
                  ivw_lme_pval = M2.2_ivw_lme$tTable[2,5],
                  egger_lm_pval = M2.2_egger_lm_pval,
                  egger_lme_pval = M2.2_egger_lme$tTable[-1,5])


##M3.2
  stage2Xboth_lmm1 = lm(y_ukb~Xhat_ukb+Xhat2_ukb)
  stage2Xboth_lm1 = summary(stage2Xboth_lmm1)
  if(length(X1snps)>1){
  mu3.2 = stage2Xboth_lmm1$fitted.values
  Gmu3.2 = t(ukb_snp_bed[,X1snps,drop=FALSE]) %*% mu3.2
  my_TEDE_Sc3.2 = my_TEDE_Sc(
                    Gmu=Gmu3.2,
                    LDcov=LDcov1,GY=GY1,YY=YY,N=n_ukb,sigY2=stage2Xboth_lm1$sigma)
  my_TEDE_aSPU3.2 = my_TEDE_aSPU(
                    Gmu=Gmu3.2,
                    LDcov=LDcov1,GY=GY1,YY=YY,N=n_ukb,sigY2=stage2Xboth_lm1$sigma)
  X1snps_mat = ukb_snp_bed[,X1snps,drop=FALSE]
  X1snps_sum_vec = apply(X1snps_mat,1,sum)
  m1.block<-list(Id = pdIdent(~X1snps_mat-1))
  stage2Xboth_egger_lme1 = tryCatch(lme(y_ukb ~ Xhat_ukb + Xhat2_ukb + X1snps_sum_vec,random = m1.block,method='ML',control = lmeControl(opt = "optim")),
                        error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
  stage2Xboth_ivw_lme1 = tryCatch(lme(y_ukb ~ Xhat_ukb + Xhat2_ukb ,random = m1.block,method='ML',control = lmeControl(opt = "optim")),
                        error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
  if(is.null(stage2Xboth_egger_lme1)){
      M3.2_egger_lme_pval=M3.2_egger_lm_pval=NULL
  }else{
  M3.2_egger_lme = summary(stage2Xboth_egger_lme1)
  M3.2_egger_lme_chisq = pchisq(t(M3.2_egger_lme$coefficients$fixed[2:3]) %*% ginv(vcov(stage2Xboth_egger_lme1)[2:3,2:3]) %*% M3.2_egger_lme$coefficients$fixed[2:3],
                            df =2 ,lower.tail=F)
  M3.2_egger_lme_pval = c(M3.2_egger_lme$tTable[-1,5],M3.2_egger_lme_chisq)

      stage2Xboth_egger_lmm1 = lm(y_ukb ~ Xhat_ukb + Xhat2_ukb + X1snps_sum_vec)
      M3.2_HC = vcovHC(stage2Xboth_egger_lmm1,type='HC0')
      M3.2_Z_pval = pchisq(stage2Xboth_egger_lmm1$coefficients[4]^2/M3.2_HC[4,4],df=1,lower.tail=F)
      M3.2_Xhat_pval = pchisq(t(stage2Xboth_egger_lmm1$coefficients[2:3]) %*% ginv(M3.2_HC[2:3,2:3]) %*% stage2Xboth_egger_lmm1$coefficients[2:3],
                                df = 2, lower.tail=F)
      M3.2_egger_lm_pval = c(M3.2_Xhat_pval,M3.2_Z_pval)
    }
  if(is.null(stage2Xboth_ivw_lme1)){
      M3.2_ivw_lme_pval = NULL
  }else{
      M3.2_ivw_lme = summary(stage2Xboth_ivw_lme1)
  M3.2_ivw_lme_chisq = pchisq(t(M3.2_ivw_lme$coefficients$fixed[2:3]) %*% ginv(vcov(stage2Xboth_ivw_lme1)[2:3,2:3]) %*% M3.2_ivw_lme$coefficients$fixed[2:3],
                            df =2 ,lower.tail=F)
  M3.2_ivw_lme_pval = c(M3.2_ivw_lme$tTable[-1,5],M3.2_ivw_lme_chisq)
  }
  }else{
   my_TEDE_Sc3.2= my_TEDE_aSPU3.2= M3.2_egger_lme_pval=M3.2_egger_lm_pval=M3.2_ivw_lme_pval=NULL
  }
  stage2Xboth_1 = list(theta=stage2Xboth_lm1$coefficients[-1,1],
                  theta_var = vcov(stage2Xboth_lm1)[-1,-1],
                  theta_pval = stage2Xboth_lm1$coefficients[-1,4],
                  sigma = stage2Xboth_lm1$sigma,
                  fstat = stage2Xboth_lm1$fstatistic,
                  aic = AIC(stage2Xboth_lmm1),
                  bic = BIC(stage2Xboth_lmm1),
                  my_TEDE_Sc=my_TEDE_Sc3.2,
                  my_TEDE_aSPU=my_TEDE_aSPU3.2,
                  ivw_lme_pval = M3.2_ivw_lme_pval,
                  egger_lm_pval = M3.2_egger_lm_pval,
                  egger_lme_pval = M3.2_egger_lme_pval)

  FINAL_RESULT = list(gene=gene_name,
                      chr=chr,
                      stage2X1 = stage2X1,
                      stage2X2_v1 = stage2X2,
                      stage2X2_v2 = stage2X2_1,
                      stage2Xboth_v1 = stage2Xboth,
                      stage2Xboth_v2 = stage2Xboth_1,
                      n2 = n_ukb,
                      p_M1 = length(X1snps),
                      p_M2.1 = length(X2snps),
                      p_M3.1 = length(X3snps),
                      fX1 = lm_stage1_X1$fstatistic,
                      fX2 = lm_stage1_X2$fstatistic,
                      AICX1 = AIC_stage1_X1,
                      BICX1 = BIC_stage1_X1,
                      r2X1 = rsq_stage1_X1,
                      adjr2X1 = adjrsq_stage1_X1,
                      AICX2 = AIC_stage1_X2,
                      BICX2 = BIC_stage1_X2,
                      r2X2 = rsq_stage1_X2,
                      adjr2X2 = adjrsq_stage1_X2
  )
  real_data_result = c(real_data_result,list(FINAL_RESULT = FINAL_RESULT))
  
}
save(real_data_result,
     file = paste("./result/LDL/TWAS_X2_bs_chr",chr,".Rdata",sep=""))

## FOR simulation
df = list(x=Xhat_ukb,x2=X2hat_ukb,y=y_ukb,
          theta1=stage2X1_lmm$coefficients[2],
          theta2=stage2X2_lmm$coefficients[2],
          theta=stage2Xboth_lmm$coefficients[-1])
save(df,file='SPDYC_LDL.Rdata')

