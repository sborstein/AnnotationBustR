library(ape)
library(seqinr)

gene.list.locs= list(
tRNA.Phe<-subset(cichlid.pos[,1:4],!is.na(cichlid.pos[,3])),
rRNA12s<-subset(cichlid.pos[,c(1:2,5:6)],!is.na(cichlid.pos[,5])),
tRNA.Val<-subset(cichlid.pos[,c(1:2,7:8)],!is.na(cichlid.pos[,7])),
rRNA16s<-subset(cichlid.pos[,c(1:2,9:10)],!is.na(cichlid.pos[,9])),
tRNA.Leu1<-subset(cichlid.pos[,c(1:2,11:12)],!is.na(cichlid.pos[,11])),
nd1<-subset(cichlid.pos[,c(1:2,13:14)],!is.na(cichlid.pos[,13])),
tRNA.Ile<-subset(cichlid.pos[,c(1:2,15:16)],!is.na(cichlid.pos[,15])),
tRNA.Gln<-subset(cichlid.pos[,c(1:2,17:18)],!is.na(cichlid.pos[,17])),
tRNA.Met<-subset(cichlid.pos[,c(1:2,19:20)],!is.na(cichlid.pos[,19])),
nd2<-subset(cichlid.pos[,c(1:2,21:22)],!is.na(cichlid.pos[,21])),
tRNA.Trp<-subset(cichlid.pos[,c(1:2,23:24)],!is.na(cichlid.pos[,23])),
tRNA.Ala<-subset(cichlid.pos[,c(1:2,25:26)],!is.na(cichlid.pos[,25])),
tRNA.Asn<-subset(cichlid.pos[,c(1:2,27:28)],!is.na(cichlid.pos[,27])),
tRNA.Cys<-subset(cichlid.pos[,c(1:2,29:30)],!is.na(cichlid.pos[,29])),
tRNA.Tyr<-subset(cichlid.pos[,c(1:2,31:32)],!is.na(cichlid.pos[,31])),
coi<-subset(cichlid.pos[,c(1:2,33:34)],!is.na(cichlid.pos[,33])),
tRNA.Ser1<-subset(cichlid.pos[,c(1:2,35:36)],!is.na(cichlid.pos[,35])),
tRNA.Asp<-subset(cichlid.pos[,c(1:2,37:38)],!is.na(cichlid.pos[,37])),
coii<-subset(cichlid.pos[,c(1:2,39:40)],!is.na(cichlid.pos[,39])),
tRNA.Lys<-subset(cichlid.pos[,c(1:2,41:42)],!is.na(cichlid.pos[,41])),
atp8<-subset(cichlid.pos[,c(1:2,43:44)],!is.na(cichlid.pos[,43])),
atp6<-subset(cichlid.pos[,c(1:2,45:46)],!is.na(cichlid.pos[,45])),
coiii<-subset(cichlid.pos[,c(1:2,47:48)],!is.na(cichlid.pos[,47])),
tRNA.Gly<-subset(cichlid.pos[,c(1:2,49:50)],!is.na(cichlid.pos[,49])),
nd3<-subset(cichlid.pos[,c(1:2,51:52)],!is.na(cichlid.pos[,51])),
tRNA.Arg<-subset(cichlid.pos[,c(1:2,53:54)],!is.na(cichlid.pos[,53])),
nd4l<-subset(cichlid.pos[,c(1:2,55:56)],!is.na(cichlid.pos[,55])),
nd4<-subset(cichlid.pos[,c(1:2,57:58)],!is.na(cichlid.pos[,57])),
tRNA.His<-subset(cichlid.pos[,c(1:2,59:60)],!is.na(cichlid.pos[,59])),
tRNA.Ser2<-subset(cichlid.pos[,c(1:2,61:62)],!is.na(cichlid.pos[,61])),
tRNA.Leu2<-subset(cichlid.pos[,c(1:2,63:64)],!is.na(cichlid.pos[,63])),
nd5<-subset(cichlid.pos[,c(1:2,65:66)],!is.na(cichlid.pos[,65])),
nd6<-subset(cichlid.pos[,c(1:2,67:68)],!is.na(cichlid.pos[,67])),
tRNA.Glu<-subset(cichlid.pos[,c(1:2,69:70)],!is.na(cichlid.pos[,69])),
cytb<-subset(cichlid.pos[,c(1:2,71:72)],!is.na(cichlid.pos[,71])),
tRNA.Thr<-subset(cichlid.pos[,c(1:2,73:74)],!is.na(cichlid.pos[,73])),
tRNA.Pro<-subset(cichlid.pos[,c(1:2,75:76)],!is.na(cichlid.pos[,75])),
dloop<-subset(cichlid.pos[,c(1:2,77:78)],!is.na(cichlid.pos[,77]))
)



for (i in sequence(length(gene.list.locs))){
  for(j in 1:dim(gene.list.locs[[i]])[1]){
    current.seq<-read.GenBank(gene.list.locs[[i]][j,2],species.names=TRUE,as.character=TRUE)#access number for read genbank
    cut.seq<-current.seq[[1]][gene.list.locs[[i]][j,3]:gene.list.locs[[i]][j,4]]#cut based on start/stop position
    seq.ids<-attr(current.seq, "species")#extract attribute-species name, from genbank takedown above
    write.fasta(sequences=cut.seq, names=seq.ids, file.out=paste(strsplit(colnames(gene.list.locs[[i]])[3],".",fixed=TRUE)[[1]][1],"fa", sep="."),open="a")#write to fasta with gene name replacing the accession# with the species name
  }
}
