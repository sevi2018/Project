setwd("/Users/emocci1/Documents/regions_1000Genomes/phasing_homozygose/LAST_may1017")

#TFB2M
start_time= Sys.time()
x<- read.table("reg16_common.ped",header=0)
n_snps=(ncol(x)-6)/2

a=seq(7,ncol(x),2)
cols_a=paste('V',a,sep='')
x$al1= apply(x[, cols_a],1,paste,collapse='')
b=seq(8,ncol(x),2)
cols_b=paste('V',b,sep='')
x$al2= apply(x[, cols_b],1,paste,collapse='')
x1=x[c(2,ncol(x)-1,ncol(x))]
colnames(x1)[1]<- "id"
x1$eq<- ifelse(x1$al1==x1$al2,1,0)###
table(x1$eq)
## 0    1
##1766  738

x1$geno<- ifelse(x1$eq==0,'het','hom')
table(x1$geno)

t<- x1

for(i in 1:n_snps){
t[,paste('snp',i,'_1',sep='')]<- substring(t$al1,i,i)
t[,paste('snp',i,'_2',sep='')]<- substring(t$al2,i,i)
i= i + 1
}


for(i in 1:n_snps){
    t[,paste('snp',i,sep='')]<- paste(t[,paste('snp',i,'_1',sep='')],t[,paste('snp',i,'_2',sep='')],sep='')
i=i+1
}

t1<- t[c(1:3,5,((ncol(t)-n_snps)+1):ncol(t))]
 
 for(i in 1:n_snps){
 t1[,paste('snp',i,'_n',sep='')]<- ifelse(t1[,paste('snp',i,sep='')]=="AA",0,
 ifelse(t1[,paste('snp',i,sep='')]=="CC",1,
 ifelse(t1[,paste('snp',i,sep='')]=="GG",2,
 ifelse(t1[,paste('snp',i,sep='')]=="TT",3,
 ifelse(t1[,paste('snp',i,sep='')]=="AC"|t1[,paste('snp',i,sep='')]=="CA",4,
 ifelse(t1[,paste('snp',i,sep='')]=="AG"|t1[,paste('snp',i,sep='')]=="GA",5,
 ifelse(t1[,paste('snp',i,sep='')]=="AT"|t1[,paste('snp',i,sep='')]=="TA",6,
 ifelse(t1[,paste('snp',i,sep='')]=="CG"|t1[,paste('snp',i,sep='')]=="GC",7,
 ifelse(t1[,paste('snp',i,sep='')]=="CT"|t1[,paste('snp',i,sep='')]=="TC",8,
 ifelse(t1[,paste('snp',i,sep='')]=="GT"|t1[,paste('snp',i,sep='')]=="TG",9,'NA'))))))))))
 }

c=seq(1:n_snps)
cols_c=paste('snp',c,'_n',sep='')
t1$all<- apply(t1[, cols_c],1,paste,collapse='')


###Select the homozygose
hom<- t1[which(t1$geno=="hom"),]
hom1=hom[c(1,2,4,ncol(hom))]
hom1$s<-1
hom2<- aggregate(hom1$s,by=list(hom1$all),FUN=sum)
names(hom2)<- c('all','nb')
hom3<- merge(hom1,hom2,by='all')
hom4=hom3[c(3,1,4,6)]
hom5=hom4[!duplicated(hom4),]
hom6<- hom5[order(hom5$nb,decreasing=TRUE),]
hom6$Var1<- LETTERS[seq(1:nrow(hom6))]
hom6$Var2<- LETTERS[seq(1:nrow(hom6))]### 7 alleles

#    al1      all geno  nb Var1 Var2
#  GCGGCGCC 21221211  hom 425    A    A
#  ATGGTAGT 03223023  hom 100    B    B
#  ACGGTAGT 01223023  hom  81    C    C
#  ACGACGCC 01201211  hom  54    D    D
#  ATGGCAGC 03221021  hom  45    E    E
#  ACGGCGGC 01221221  hom  26    F    F
#  GCGGCGGC 21221221  hom   7    G    G
###Create all possible combination of heterozygouse given the known homoz alleles

cc<- expand.grid(hom6$Var1,hom6$Var2)
cc$eq<- ifelse(cc$Var1==cc$Var2,0,1)
cc1<- cc[which(cc$eq==1),]
cc1$eq<- seq(1:nrow(cc1))

hz1=hom6[c(5,1,2)]
hz2=hom6[c(6,1,2)]
hz1b= merge(hz1,cc1,by='Var1')
colnames(hz2)[2]='al2'
hz3=merge(hz1b,hz2,by='Var2')
hz4=hz3[c(2,1,3,4,6,7,5)]
hz5=hz4[order(hz4$eq),]

for(i in 1:n_snps){
    hz5[,paste('snp',i,'_1',sep='')]<- substring(hz5$all.x,i,i)
    hz5[,paste('snp',i,'_2',sep='')]<- substring(hz5$all.y,i,i)
    i= i + 1
}

for(i in 1:n_snps){
    hz5[,paste('snp',i,'g',sep='')]<- ifelse(hz5[,paste('snp',i,'_1',sep='')]==0 & hz5[,paste('snp',i,'_2',sep='')]==0,0,
    ifelse(hz5[,paste('snp',i,'_1',sep='')]==1 & hz5[,paste('snp',i,'_2',sep='')]==1,1,
    ifelse(hz5[,paste('snp',i,'_1',sep='')]==2 & hz5[,paste('snp',i,'_2',sep='')]==2,2,
    ifelse(hz5[,paste('snp',i,'_1',sep='')]==3 & hz5[,paste('snp',i,'_2',sep='')]==3,3,
    ifelse(hz5[,paste('snp',i,'_1',sep='')]==0 & hz5[,paste('snp',i,'_2',sep='')]==1|hz5[,paste('snp',i,'_1',sep='')]==1 & hz5[,paste('snp',i,'_2',sep='')]==0,4,
    ifelse(hz5[,paste('snp',i,'_1',sep='')]==0 & hz5[,paste('snp',i,'_2',sep='')]==2|hz5[,paste('snp',i,'_1',sep='')]==2 & hz5[,paste('snp',i,'_2',sep='')]==0,5,
    ifelse(hz5[,paste('snp',i,'_1',sep='')]==0 & hz5[,paste('snp',i,'_2',sep='')]==3|hz5[,paste('snp',i,'_1',sep='')]==3 & hz5[,paste('snp',i,'_2',sep='')]==0,6,
    ifelse(hz5[,paste('snp',i,'_1',sep='')]==1 & hz5[,paste('snp',i,'_2',sep='')]==2|hz5[,paste('snp',i,'_1',sep='')]==2 & hz5[,paste('snp',i,'_2',sep='')]==1,7,
    ifelse(hz5[,paste('snp',i,'_1',sep='')]==1 & hz5[,paste('snp',i,'_2',sep='')]==3|hz5[,paste('snp',i,'_1',sep='')]==3 & hz5[,paste('snp',i,'_2',sep='')]==1,8,
    ifelse(hz5[,paste('snp',i,'_1',sep='')]==2 & hz5[,paste('snp',i,'_2',sep='')]==3|hz5[,paste('snp',i,'_1',sep='')]==3 & hz5[,paste('snp',i,'_2',sep='')]==2,9,'NA'))))))))))
    i=i+1
}

d=seq(1:n_snps)
cols_d=paste('snp',c,'g',sep='')
hz5$all<- apply(hz5[, cols_d],1,paste,collapse='')

hz6=hz5[c(1,2,ncol(hz5))]
hz6_= hz6[!duplicated(hz6$all),]
hz6_$geno='het'
homz= hom6[c(5,6,2,3)]
als= rbind(homz,hz6_)
table(als$geno)

# het hom
# 21  7
alleles=hom6[c(1,5)]
al1=merge(als,alleles,by='Var1')
names(alleles)=c('al2','Var2')
al2=merge(al1,alleles,by='Var2')
al3=al2[c(2,1,5,6,3,4)]
al3$g=paste(al3$Var1,al3$Var2,sep='')

### MATCH THE HET IN SILICO BUILT FROM HOMOZYGOUSE ALLELES WITH THE CODE OF THE SAMPLES
t2= t1[c(1:4,ncol(t1))]
t2$s= seq(1:nrow(t2))
t3= merge(t2,al3,by='all',all=TRUE)
t4=t3[order(t3$s),]
t5= t4[which(t4$id!='NA'),]
####
#######IDs identified by homozygous alleles
id0=t5[c(2,1,ncol(t5))]
id0_=id0[which(id0$g!='NA'),]
id0_$id_hom=1
###How many genotype am I missing?
unk= t5[which(is.na(t5$g=='NA')),]##124 (4.9%)
bo5=unk


##FIND NEW ALLELES::LOOP1
sco=bo5[c(1)]
sco$x=1
sc1=aggregate(sco$x,by=list(sco$all),FUN=sum)
sc2=sc1[order(-sc1$x),]
sc3=sc2[which(sc2$x>1),]

n1=hom6[c(1,2)]

###Decompose the code into the possible alleles starting from the most common one
x=1:nrow(sc3)
for (val in x) {
    one=as.data.frame(sc3[val,])

for(i in 1:n_snps){
    one[,paste('snp',i,sep='')]<- substring(one$Group.1,i,i)
    i= i + 1
}
sapply(one,class)##"character"
cols.num=c(3:ncol(one))
one[cols.num]=sapply(one[cols.num],as.numeric)
sapply(one,class)##"numeric"
o=one[,3:ncol(one)]
o1=as.data.frame(lapply(o,function(x) x>=4))
o2=as.data.frame(lapply(o,function(x) x<=3))
all=rbind(o,o1,o2)
al=as.data.frame(t(all))
names(al)=c('als','het','hom')
het=al[which(al$het==1),]

het$a1=ifelse(het$als==4,0,
ifelse(het$als==5,0,
ifelse(het$als==6,0,
ifelse(het$als==7,1,
ifelse(het$als==8,1,
ifelse(het$als==9,2,'NA'))))))

het$a2=ifelse(het$als==4,1,
ifelse(het$als==5,2,
ifelse(het$als==6,3,
ifelse(het$als==7,2,
ifelse(het$als==8,3,
ifelse(het$als==9,3,'NA'))))))

l=as.data.frame(t(het))
l1=l[4:5,]
het1=expand.grid(l1)
names(het1)=names(l)

hom=al[which(al$hom==1),]
hm=as.data.frame(t(hom))
hm1=hm[1,]
hm2=hm1[rep(seq_len(nrow(hm1)),each=nrow(het1)),]
hm3=cbind(hm2,het1)
hm4=hm3[,order(names(hm3))]
cols_e=paste('snp',c,sep='')
hm4$new_als<- apply(hm4[, cols_e],1,paste,collapse='')
hm5=hm4[!duplicated(hm4$new_als),]
new=hm5[c(ncol(hm5))]
names(new)='all'
new$x='new'
hm6=merge(hom6,new,by='all')

if (nrow(hm6) > 1 |nrow(hm6) == 0){
    next
}

## subctract the combined code to the known allele and select the new allele
old=hm6[c(1,2)]
for(i in 1:n_snps){
old[,paste('snp',i,sep='')]<- substring(old$all,i,i)
i= i + 1
}
cb=one[c(1,3:ncol(one))]
cols_f=paste('s',c,sep='')
names(cb)[2:(n_snps+1)]<- cols_f
names(cb)[1]<- 'code'
oo=cbind(cb,old)

for(i in 1:n_snps){
    oo[,paste('snp',i,'n',sep='')]<- ifelse(oo[,paste('s',i,sep='')]==0 & oo[,paste('snp',i,sep='')]==0,0,
    ifelse(oo[,paste('s',i,sep='')]==1 & oo[,paste('snp',i,sep='')]==1,1,
    ifelse(oo[,paste('s',i,sep='')]==2 & oo[,paste('snp',i,sep='')]==2,2,
    ifelse(oo[,paste('s',i,sep='')]==3 & oo[,paste('snp',i,sep='')]==3,3,
    ifelse(oo[,paste('s',i,sep='')]==4 & oo[,paste('snp',i,sep='')]==1,0,
    ifelse(oo[,paste('s',i,sep='')]==4 & oo[,paste('snp',i,sep='')]==0,1,
    ifelse(oo[,paste('s',i,sep='')]==5 & oo[,paste('snp',i,sep='')]==2,0,
    ifelse(oo[,paste('s',i,sep='')]==5 & oo[,paste('snp',i,sep='')]==0,2,
    ifelse(oo[,paste('s',i,sep='')]==6 & oo[,paste('snp',i,sep='')]==3,0,
    ifelse(oo[,paste('s',i,sep='')]==6 & oo[,paste('snp',i,sep='')]==0,3,
    ifelse(oo[,paste('s',i,sep='')]==7 & oo[,paste('snp',i,sep='')]==1,2,
    ifelse(oo[,paste('s',i,sep='')]==7 & oo[,paste('snp',i,sep='')]==2,1,
    ifelse(oo[,paste('s',i,sep='')]==8 & oo[,paste('snp',i,sep='')]==1,3,
    ifelse(oo[,paste('s',i,sep='')]==8 & oo[,paste('snp',i,sep='')]==3,1,
    ifelse(oo[,paste('s',i,sep='')]==9 & oo[,paste('snp',i,sep='')]==2,3,
    ifelse(oo[,paste('s',i,sep='')]==9 & oo[,paste('snp',i,sep='')]==3,2,9))))))))))))))))
    i=i+1
}

cols_g=paste('snp',c,'n',sep='')
oo$new_als<- apply(oo[, cols_g],1,paste,collapse='')
###Build the new allele from the code
for(i in 1:n_snps){
oo[,paste('known',i,sep='')]<- substring(oo$al1,i,i)
i= i + 1
}
j=ncol(oo)-n_snps
z=(ncol(oo)-n_snps)-(n_snps*2)
o1=oo[c((j+1):ncol(oo),z:(j-1))]

for(i in 1:n_snps){
    o1[,paste('new',i,sep='')]<- ifelse(o1[,paste('snp',i,sep='')]==0 & o1[,paste('snp',i,'n',sep='')]==0 & o1[,paste('known',i,sep='')]=='A','A',
    o1[,paste('new',i,sep='')]<- ifelse(o1[,paste('snp',i,sep='')]==0 & o1[,paste('snp',i,'n',sep='')]==1 & o1[,paste('known',i,sep='')]=='A','C',
    o1[,paste('new',i,sep='')]<- ifelse(o1[,paste('snp',i,sep='')]==1 & o1[,paste('snp',i,'n',sep='')]==0 & o1[,paste('known',i,sep='')]=='C','A',
    o1[,paste('new',i,sep='')]<- ifelse(o1[,paste('snp',i,sep='')]==0 & o1[,paste('snp',i,'n',sep='')]==2 & o1[,paste('known',i,sep='')]=='A','G',
    o1[,paste('new',i,sep='')]<- ifelse(o1[,paste('snp',i,sep='')]==2 & o1[,paste('snp',i,'n',sep='')]==0 & o1[,paste('known',i,sep='')]=='G','A',
    o1[,paste('new',i,sep='')]<- ifelse(o1[,paste('snp',i,sep='')]==0 & o1[,paste('snp',i,'n',sep='')]==3 & o1[,paste('known',i,sep='')]=='A','T',
    o1[,paste('new',i,sep='')]<- ifelse(o1[,paste('snp',i,sep='')]==3 & o1[,paste('snp',i,'n',sep='')]==0 & o1[,paste('known',i,sep='')]=='T','A',
    o1[,paste('new',i,sep='')]<- ifelse(o1[,paste('snp',i,sep='')]==1 & o1[,paste('snp',i,'n',sep='')]==1 & o1[,paste('known',i,sep='')]=='C','C',
    o1[,paste('new',i,sep='')]<- ifelse(o1[,paste('snp',i,sep='')]==1 & o1[,paste('snp',i,'n',sep='')]==2 & o1[,paste('known',i,sep='')]=='C','G',
    o1[,paste('new',i,sep='')]<- ifelse(o1[,paste('snp',i,sep='')]==2 & o1[,paste('snp',i,'n',sep='')]==1 & o1[,paste('known',i,sep='')]=='G','C',
    o1[,paste('new',i,sep='')]<- ifelse(o1[,paste('snp',i,sep='')]==1 & o1[,paste('snp',i,'n',sep='')]==3 & o1[,paste('known',i,sep='')]=='C','T',
    o1[,paste('new',i,sep='')]<- ifelse(o1[,paste('snp',i,sep='')]==3 & o1[,paste('snp',i,'n',sep='')]==1 & o1[,paste('known',i,sep='')]=='T','C',
    o1[,paste('new',i,sep='')]<- ifelse(o1[,paste('snp',i,sep='')]==2 & o1[,paste('snp',i,'n',sep='')]==2 & o1[,paste('known',i,sep='')]=='G','G',
    o1[,paste('new',i,sep='')]<- ifelse(o1[,paste('snp',i,sep='')]==2 & o1[,paste('snp',i,'n',sep='')]==3 & o1[,paste('known',i,sep='')]=='G','T',
    o1[,paste('new',i,sep='')]<- ifelse(o1[,paste('snp',i,sep='')]==3 & o1[,paste('snp',i,'n',sep='')]==2 & o1[,paste('known',i,sep='')]=='T','G',
    o1[,paste('new',i,sep='')]<- ifelse(o1[,paste('snp',i,sep='')]==3 & o1[,paste('snp',i,'n',sep='')]==3 & o1[,paste('known',i,sep='')]=='T','T','NA'))))))))))))))))
}
cols_h=paste('new',c,sep='')
o1$al1=apply(o1[, cols_h],1,paste,collapse='')
cols_i=paste('snp',c,'n',sep='')
o1$all=apply(o1[, cols_i],1,paste,collapse='')




######## Put the new allele together with the ones identified from homoZ and repeat the identification process
nh=n1[c(1,2)]
new=o1[c((ncol(o1)-1),ncol(o1))]
n1=rbind(nh,new)
n1=n1[!duplicated(n1),]
LETTERS100=c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))
n1$Var1<- LETTERS100[seq(1:nrow(n1))]
n1$Var2<- LETTERS100[seq(1:nrow(n1))]

cc<- expand.grid(n1$Var1,n1$Var2)
cc$eq<- ifelse(cc$Var1==cc$Var2,0,1)
cc1<- cc[which(cc$eq==1),]
cc1$eq<- seq(1:nrow(cc1))

n2= merge(n1,cc1,by='Var1')
n3=n2[c(3,1,5,6)]
colnames(n3)[3]='Var2'
n4=merge(n3,n1,by='Var2')
n5=n4[order(n4$eq),]
n6=n5[c(4,1,6,3,2)]
names(n6)=c('Nb','Var1','al1','Var2','al2')

###Build new combination of alleles using the new allele

for(i in 1:n_snps){
    n6[,paste('snp',i,'_1',sep='')]<- substring(n6$al1,i,i)
    n6[,paste('snp',i,'_2',sep='')]<- substring(n6$al2,i,i)
    i= i + 1
}

for(i in 1:n_snps){
    n6[,paste('snp',i,'g',sep='')]<- ifelse(n6[,paste('snp',i,'_1',sep='')]==0 & n6[,paste('snp',i,'_2',sep='')]==0,0,
    ifelse(n6[,paste('snp',i,'_1',sep='')]==1 & n6[,paste('snp',i,'_2',sep='')]==1,1,
    ifelse(n6[,paste('snp',i,'_1',sep='')]==2 & n6[,paste('snp',i,'_2',sep='')]==2,2,
    ifelse(n6[,paste('snp',i,'_1',sep='')]==3 & n6[,paste('snp',i,'_2',sep='')]==3,3,
    ifelse(n6[,paste('snp',i,'_1',sep='')]==0 & n6[,paste('snp',i,'_2',sep='')]==1|n6[,paste('snp',i,'_1',sep='')]==1 & n6[,paste('snp',i,'_2',sep='')]==0,4,
    ifelse(n6[,paste('snp',i,'_1',sep='')]==0 & n6[,paste('snp',i,'_2',sep='')]==2|n6[,paste('snp',i,'_1',sep='')]==2 & n6[,paste('snp',i,'_2',sep='')]==0,5,
    ifelse(n6[,paste('snp',i,'_1',sep='')]==0 & n6[,paste('snp',i,'_2',sep='')]==3|n6[,paste('snp',i,'_1',sep='')]==3 & n6[,paste('snp',i,'_2',sep='')]==0,6,
    ifelse(n6[,paste('snp',i,'_1',sep='')]==1 & n6[,paste('snp',i,'_2',sep='')]==2|n6[,paste('snp',i,'_1',sep='')]==2 & n6[,paste('snp',i,'_2',sep='')]==1,7,
    ifelse(n6[,paste('snp',i,'_1',sep='')]==1 & n6[,paste('snp',i,'_2',sep='')]==3|n6[,paste('snp',i,'_1',sep='')]==3 & n6[,paste('snp',i,'_2',sep='')]==1,8,
    ifelse(n6[,paste('snp',i,'_1',sep='')]==2 & n6[,paste('snp',i,'_2',sep='')]==3|n6[,paste('snp',i,'_1',sep='')]==3 & n6[,paste('snp',i,'_2',sep='')]==2,9,'NA'))))))))))
    i=i+1
}

cols_l=paste('snp',c,'g',sep='')
n6$all=apply(n6[, cols_l],1,paste,collapse='')

n7=n6[c(2,4,ncol(n6))]
n7_= n7[!duplicated(n7$all),]
n7_$geno='het'
homz= n1[c(2,3,4)]
homz$geno='hom'
als= rbind(homz,n7_)
table(als$geno)

# het hom
# 28  8

als$g=paste(als$Var1,als$Var2,sep='')

###Try to identify the alleles of people with still no genotype idenityfied
bo1=bo5[c(2,1,6)]
names(bo1)=c('id','all','s')
bo1$s= seq(1:nrow(bo1))
bo2= merge(bo1,als,by='all',all=TRUE)
bo3=bo2[order(bo2$s),]
bo4= bo3[which(bo3$id!='NA'),]
bo4$id_het=1
bo4_=bo4[which(bo4$g!='NA'),]
bo6=bo4_[c(2,1,7,8)]

bo5= bo4[which(is.na(bo4$g=='NA')),]##66
bo5_=bo5[c(2,1,7,8)]

####LOOP1:GCGGTGGC 21223221    H
write.table(bo6,paste('id',val,'.txt', sep=''),row.names=F, col.names=T, quote=F)
}
write.table(id0_,'id0',row.names=F,col.names=T,quote=F)
write.table(bo5_,'ids_not_phased.txt',row.names=F,col.names=T,quote=F)
end_time <- Sys.time()
end_time - start_time



##TERMINAL
for i in $(seq 1 54); do sed '1d' id"$i".txt|cat >>SORCS2_CHR4; done
wc -l SORCS2_CHR4
cat id0 SORCS2_CHR4>>chr4_SORCS2.txt
wc -l chr4_SORCS2.txt
wc -l ids_not_phased.txt
for i in $(seq 1 22); do rm id"$i".txt; done
rm id0
