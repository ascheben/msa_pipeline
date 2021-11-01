args<-commandArgs(TRUE)
GERP_clade1<-read.table(args[1])
GERP_clade2<-read.table(args[2])
coorderate.position<-read.table(args[3],header =T)

GERP_clade1$genome.pos<-paste(GERP_clade1[,1],GERP_clade1[,2],sep="_")
GERP_clade2$genome.pos=paste(GERP_clade2[,1],GERP_clade2[,2],sep="_")

coorderate.position$Target_genome.pos=paste(coorderate.position$Target_chr,coorderate.position$Target_position,sep="_")
coorderate.position$Quary_genome.pos=paste(coorderate.position$Quary_chr,coorderate.position$Quary_position,sep="_")

index_clade1=match(coorderate.position$Target_genome.pos,GERP_clade1$genome.pos)
coorderate.position$GERP_clade1.GERP=GERP_clade1[index_clade1,4]

index_clade2=match(coorderate.position$Quary_genome.pos,GERP_clade2$genome.pos)
coorderate.position$GERP_clade2.GERP=GERP_clade2[index_clade2,4]


coorderate.position_final=coorderate.position[(!is.na(coorderate.position$GERP_clade1.GERP) & !is.na(coorderate.position$GERP_clade2.GERP)),]
shared.GERP.loci=sum(!is.na(coorderate.position_final$GERP_clade1.GERP) & !is.na(coorderate.position_final$GERP_clade2.GERP) )


shared.positive.GERP.loci=sum((coorderate.position_final$GERP_clade1.GERP>0 & coorderate.position_final$GERP_clade2.GERP>0),na.rm=T)

median.shared.GERP1.loci=median(coorderate.position_final$GERP_clade1.GERP[(!is.na(coorderate.position_final$GERP_clade1.GERP) & !is.na(coorderate.position_final$GERP_clade2.GERP))])
median.shared.GERP2.loci=median(coorderate.position_final$GERP_clade1.GERP[(!is.na(coorderate.position_final$GERP_clade1.GERP) & !is.na(coorderate.position_final$GERP_clade2.GERP))])

GERP_cor.r=cor(coorderate.position_final$GERP_clade1.GERP,coorderate.position_final$GERP_clade2.GERP,use="pairwise.complete.obs")

Clade1.GERP.loci=nrow(GERP_clade1)
Clade2.GERP.loci=nrow(GERP_clade2)
two.clade.gerp.estimate=cbind(Clade1.GERP.loci,Clade2.GERP.loci,shared.GERP.loci,GERP_cor.r,shared.positive.GERP.loci)

write.table(two.clade.gerp.estimate,args[4],quote=F,sep="\t",row.name=F)
