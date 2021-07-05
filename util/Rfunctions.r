





#################### functions for calculate slopes

calculate_piNpis=function(outdir,rdir=NULL,filter=T,numChr=NULL,miss.cutoff=0.5,allSite=F,resample=T, common.list=NULL, eval=1e-20,bscore=200,qcov=90,qlen=30, base=NULL,bad_snp.cutoff=10,write.genes=NULL, size=NULL){
	data.process<-process_data(outdir, rdir,filter, numChr,miss.cutoff,allSite,resample, common.list, eval, bscore,qcov, qlen, bad_snp.cutoff)
	new.data<-data.process[[1]]
	if(!is.null(size)){
		if(size>nrow(new.data)){
			stop("not enough genes to subsample")
		}else{
			new.data<-new.data[sample(1:nrow(new.data), size=size, rep=F),]
		}
	}
	new.data$Pi0<-ifelse(new.data$L0.correct==0, NA, new.data$Pi0/new.data$L0.correct) 
	new.data$Pi4<-ifelse(new.data$L4.correct==0, NA, new.data$Pi4/new.data$L4.correct)
	
	
	if(is.null(numChr)){numChr<-data.process[[2]]}
	
	num.genes<-nrow(new.data)
	num.total_variants<-sum(new.data$Num_SNP)
	if(allSite){
		average_miss_rates<-sum(new.data$bad_snp)/sum(new.data$Length)*100
	}else{
		average_miss_rates<-sum(new.data$bad_snp)/sum(new.data$Num_SNP)*100
	}
	
	num.good_variants<-sum(new.data$f0+new.data$f4+new.data$others) ## f0 == snp0 f4==snp4
	num.f0.variants<-sum(new.data$SNP0) 
	num.f4.variants<-sum(new.data$SNP4)
	new.data= new.data[,c(1:9,11,12,14:23)]
	x<-summary(new.data)
	res<-as.numeric(gsub('\\s', '', unlist(strsplit(x[4,c(3:5,8,9)], ':'))[seq(2,10,2)]))
	res<-c(num.total_variants,average_miss_rates, num.good_variants, num.genes, res[1:3], num.f0.variants,num.f4.variants, res[4:5])
	names(res)<-c("allSNP","Miss","usedSNP", "genes","L","L0","L4","SNP0","SNP4","Pi0","Pi4")
	print(res)
	print("Summmary Table")
	print(x)
	### create SFS and DFE input	f
	freq4<-read.table(paste(outdir,'freq_4.txt',sep='/'),header=T)
	freq0<-read.table(paste(outdir,'freq_0.txt',sep='/'),header=T)
	freq4<-freq4[freq4$Transcript%in%new.data$Transcript,]
	freq0<-freq0[freq0$Transcript%in%new.data$Transcript,]
	sfs4=summary(factor(freq4$Freq))
	sfs0=summary(factor(freq0$Freq))

	if(length(sfs0)< numChr-1){
		sfs0 <- add_null(sfs0, numChr)
	}
	if(length(sfs4)< numChr-1){
		sfs4 <- add_null(sfs4,numChr)
	}
	sfs<-rbind(sfs4, sfs0)
#	sfs<-rbind(summary(factor(freq4$Freq)), summary(factor(freq0$Freq)))
	sfs<-as.matrix(sfs)
	pdf(file='sfs.pdf')
	barplot(sfs,beside=T,args.legend=list('topright'), legend.text = c('4-fold','0-fold'))	
	dev.off()
	if(is.null(base)){
		comment=paste('data',res[4],'genes',sep='_')
		dfe.file<-paste(outdir, 'dfe_input.txt',sep='/')
	}else{
		comment=paste(base, res[4],'genes',sep='_')
		dfe.file<-paste(outdir,'/', base,'_dfe_input.txt',sep='')
	}
	cat(paste(1,1,numChr,sep="\t"),"\n",sep='',file=dfe.file)
	cat(sfs[1,],sum(new.data$L4.correct),sep="\t",append=T, file=dfe.file)
	cat("\n", sep='',append=T,file=dfe.file)
	cat(sfs[2,],sum(new.data$L0.correct),sep="\t",append=T, file=dfe.file)
	cat("\n", sep='',append=T,file=dfe.file)
	print("SFS")
	print(sfs)
	if(is.null(write.genes)==F){
		write.table(new.data,file=paste(outdir,write.genes,sep='/'),quote=F,row.names=F,col.names=T,sep="\t")
	}
	return(res)
}

add_null=function(sfs, n){
	step = 1/n
	bins = seq(step,1-step,step)
	miss_freq = bins[!bins%in%names(sfs)]
	miss = numeric(length(miss_freq))
	miss=as.data.frame(t(miss))
	names(miss)=miss_freq
	sfs<-cbind(t(sfs),miss)
	index= order(names(sfs))
	sfs=sfs[,index]
	return(sfs)
}

process_data=function(outdir,rdir=NULL,filter=T,numChr=NULL,miss.cutoff=0.5,allSite=F,resample=T, common.list=NULL, eval=1e-20,bscore=200,qcov=90, qlen=30, bad_snp.cutoff=10){
	setwd(outdir)
	### get files
	files=list.files(path=outdir, pattern='txt')
	if(is.null(numChr)){
		m<-gregexpr("Ind\\d+", files,ignore.case = T)
		numChr<-as.integer(sub('Ind', '', unlist(regmatches(files,m))[[1]],ignore.case = T))
	}
	freq<-grep('frequency_DAF', files, value=T)
	if(length(freq)==0){
		freq<-grep('frequency.txt', files, value=T)
		warning(paste("no DAF freq file found, going to use the MAF freq file"))
	}
	
	data_file<-grep('Rout.txt$', grep('log',grep('piNpiS', files,value=T, ignore.case=T),value=T,invert=T, ignore.case=T), value=T,invert=T, ignore.case=T)
	data<-read.table(data_file, header=T)
	if(grepl('diploid', data_file, ignore.case=T)){
		numChr=numChr*2
	}
	miss<-read.table(grep('missing_ratio', files, value=T, ignore.case=T),header=T)
	snp_file<-grep('num_snps', files, ignore.case=T,value=T)

	snp<-read.table(snp_file,header=T)
	miss<-merge(miss,snp,by="Transcript")
	data<-merge(data,miss,by="Transcript")
	data$bad_snp<-data$Num_N # data$Num_N is failed sequencing sites, data$fN is the sites has ambiguous fold
	
	if(allSite){
		data$L4.correct<-data$L4
		data$L0.correct<-data$L0
		print('Data is treated as full sites')
	}else{
		data$L4.correct<-ifelse(data$Num_SNP==0 | data$Num_SNP==data$Num_N, data$L4, (data$L4-data$SNP4)*(data$Num_SNP-data$Num_N)/data$Num_SNP+data$SNP4)
		data$L0.correct<-ifelse(data$Num_SNP==0 | data$Num_SNP==data$Num_N, data$L0, (data$L0-data$SNP0)*(data$Num_SNP-data$Num_N)/data$Num_SNP+data$SNP0)
	}
	# data$Pi0<-ifelse(data$L0.correct==0, NA, data$Pi0/data$L0.correct)
	# data$Pi4<-ifelse(data$L4.correct==0, NA, data$Pi4/data$L4.correct)
	data$Gene_TD0=NULL
	data$Gene_TD4=NULL
	for(i in 1:nrow(data)){
		data$Gene_TD0[i]=add_tajimasD_per_gene(data$Pi0[i], data$SNP0[i], numChr)
		data$Gene_TD4[i]=add_tajimasD_per_gene(data$Pi4[i], data$SNP4[i], numChr)
	}
		
	if(!is.null(common.list)){
		print('Selected results based on the gene list')
		geneList<-read.table(common.list,header=T)
		data<-data[data$Transcript%in%geneList$Transcript,]
	}else if(filter){
		print('Performing filtering using blast results, and only keep the longest transcripts')
		if(is.null(rdir)){stop("require reference directory")}else{
			r.files=list.files(path=rdir, pattern='txt',full.names=T)
			try(blastp<-read.table(grep('blastp', r.files, value=T, ignore.case=T)), silent=T)
			try(part<-read.table(grep('partial', r.files,value=T, ignore.case=T),sep="\t"),silent=T)
			try(longest<-read.table(grep('longest', r.files, value=T,ignore.case=T),header=T),silent=T)
			
			if(exists('blastp')){
				best.hit<-blastp[blastp$V5<=eval & blastp$V6>=bscore & blastp$V2>qlen,]
				data<-data[data$Transcript%in%best.hit$V1,]
			}else{
				warning("No blastp file exists")
			}
			if(exists('part')){
				data<-data[!data$Transcript%in%part$V1,]
			}else{
				warning("no partial CDS info found")
			}
			if(exists('longest')){
				data<-data[data$Transcript%in%longest$Transcript,]
			}else{
				warning("no longest CDS info found")
			}
			data<-data[!is.na(data$Pi0) & !is.na(data$Pi4),] # remove those bad genes
			data$used_snp<-data$SNP0+data$SNP4
			
			miss_rate<-ifelse(data$used_snp==0, 0, data$bad_snp/data$Num_SNP)
			#if(allSite){miss_rate=data$bad_snp/data$Length} 		 ### if all sites are reported then we use gene length
			if(allSite){
				filter.condition=data$bad_snp<=bad_snp.cutoff & data$bad_snp/data$Length <= miss.cutoff #### <---- filtering critering for bad snp
			}else{
				filter.condition=data$bad_snp<=bad_snp.cutoff & miss_rate<=miss.cutoff 
			}
			
			if(resample){
				print('Performing resampling based on original of non-polymorphic/polymorphic ratio')
				# if(allSite){
				# 	filter.condition=data$bad_snp<=bad_snp.cutoff | data$used_snp>=data$bad_snp
				# }else{
				# 	filter.condition=miss_rate<=miss.cutoff
				# }
				
				new.data.poly<-data[data$used_snp>0 & filter.condition ,]
				num.null.snp.genes<-round(sum(data$used_snp==0)*sum(data$used_snp[filter.condition]>0)/ sum(data$used_snp>0),dig=0)
				
				if(num.null.snp.genes <= sum(filter.condition & data$used_snp==0)){
					new.data.null<-data[sample(which(filter.condition & data$used_snp==0), num.null.snp.genes, rep=F),]
				}else{
					new.data.null<-data[filter.condition & data$used_snp==0, ]
					num.poly.genes<- round(sum(data$used_snp>0)/sum(data$used_snp==0) *  nrow(new.data.null),dig=0)
					new.data.poly<-new.data.poly[sample(1:nrow(new.data.poly), num.poly.genes,rep=F),]
				}
				new.data<-rbind(new.data.poly, new.data.null)
				data<-new.data						
			}else{
				# if(allSite){
				# 	data<-data[data$bad_snp<=bad_snp.cutoff,]
				# }else{
				# 	data<-data[miss_rate<=miss.cutoff,]
				# }
				data<-data[filter.condition,]
			}
		}		
	}
	return(list(data,numChr))
}


cal_piNpiS_by_bins_numSNPs=function(outdir,rdir=NULL,filter=T,numChr=NULL,common.list=NULL,num.bins=10,miss.cutoff=0.5,base=NULL,dfe=F,allSite=F,resample=T,rec.file=NULL,sort.key='pi',GC=NULL, eval=1e-20,bscore=200,qcov=90,qlen=30, bad_snp.cutoff=10){
	data.process<-process_data(outdir, rdir,filter, numChr,miss.cutoff,allSite,resample, common.list, eval,bscore, qcov,qlen, bad_snp.cutoff)
	data<-data.process[[1]]

	if(is.null(numChr)){numChr<-data.process[[2]]}

	# if(!is.null(common.list)){
	# 	geneList<-read.table(common.list,header=T)
	# 	data<-data[data$Transcript%in%geneList$Transcript,]
	# }

	# data$D<-add_tajimasD(data$Pi0+data$Pi4, data$SNP0+data$SNP4, data$L0.correct+data$L4.correct,numChr) ## <-- TD for syn + nonsyn
	
	# data$D<-add_tajimasD(data$Pi0, data$SNP0, data$L0.correct,numChr) ## <-- TD for nonsyn
	
	freq4<-read.table(paste(outdir,'freq_4.txt',sep='/'),header=T)
	freq4<-freq4[freq4$Transcript%in%data$Transcript,]
	freq0<-read.table(paste(outdir,'freq_0.txt',sep='/'),header=T)
	freq0<-freq0[freq0$Transcript%in%data$Transcript,]
	num.p4<-sum(data$SNP4)
	num.L4<-sum(data$L4.correct)
	
	### Fay and Wu's H
	# data$FWH<-add_fayWu_H_norm(data$Pi0+data$Pi4,  data$SNP0+data$SNP4, headdata$L0.correct+data$L4.correct,data$Transcript, numChr, rbind(freq4, freq0))
	
	## get number of snp for each group
	num.groups<-get_num_pi4_site_hyper(num.p4, num.L4)

	g1.index<-sample(1: num.p4, num.groups[1],rep=F)
	g2.index<-sample(seq(1,num.p4)[-g1.index], num.groups[2],rep=F)
	g3.index<-seq(1,num.p4)[-c(g1.index,g2.index)]

	g1<-freq4[g1.index,]
	g2<-freq4[g2.index,]
	g3<-freq4[g3.index,]
	
	data<-add_thetaW_by_group(data,g1,numChr,1)
	# data<-add_pi_by_group(data,g1,numChr,1)
	data<-add_pi_by_group(data,g2,numChr,2)
	data<-add_pi_by_group(data,g3,numChr,3)
	# data<-add_thetaW_by_group(data,g2,numChr,2)
	# data<-add_thetaW_by_group(data,g3,numChr,3)
	data$Pi_g1<-data$Pi_g1/(data$L4.correct/3)
	data<-data[!is.na(data$Pi_g1),] ## sometimes it creates NA

	### calculate number of syn SNPs for each g1 categories
	# num.snp4.p1<-summary(g1$Transcript,maxsum=1e7)
	num.snp4.p1<-summary(freq4$Transcript,maxsum=1e7)
	numSNP4<-data.frame(Transcript=names(num.snp4.p1), NumSNP4=num.snp4.p1)
	data$NumSNP4<-0
	numSNP4<-rbind(numSNP4, data[,c("Transcript","NumSNP4")])
	numSNP4<-tapply(numSNP4$NumSNP4,numSNP4$Transcript,sum)
	numSNP4<-data.frame(Transcript=names(numSNP4), NumSNP4=numSNP4)
	data<-data[,-ncol(data)]
	data<-merge(data, numSNP4,by='Transcript')

	## ranking genes
	data<-data[sample(1:nrow(data),size=nrow(data),rep=F),] # the add_pi_by_group function sort the pi values in some way, so we need to reshuffle the values
	
	if(!is.null(GC)){
		GC <-read.table(GC,header=T)
		data<-merge(data,GC,by='Transcript')
	}
	if(!is.null(rec.file)){
		gff<-read.table(rec.file,header=T)
		data<-merge(data,gff,by='Transcript')
	}
	if(sort.key == 'rec'){
		if(is.null(rec.file)) stop('No rec file specified')
		data<-data[order(data$rec),]
	}else if(sort.key =='GC'){
		if(is.null(GC)) stop('No GC file specified')
		data<-data[order(data$GC3),]
	}else{
		data<-data[order(data$Pi_g1),]
	}

	
	## rank by number of genes
	print(paste(nrow(data),"genes were investigated in total after filtering"))

	data$rank=0
	step<-round(seq(1,sum(data$NumSNP4),length=num.bins+1))[-1]
	index.bin=1
	sum.snp.bin=0
	for(i in 1:nrow(data)){
		sum.snp.bin=sum.snp.bin+data$NumSNP4[i]
		if(sum.snp.bin<=step[index.bin]){
			data$rank[i]=index.bin
		}else{
			index.bin=index.bin+1
			data$rank[i]=index.bin
		}
	}
	wt<-summary(factor(data$rank))
	pi0<-tapply(data$Pi0,data$rank,sum)/tapply(data$L0.correct, data$rank,sum)
	pi4.g2<-tapply(data$Pi_g2,data$rank,sum)/tapply(data$L4.correct,data$rank,sum)*3
	pi0.pi4<-pi0/pi4.g2
	pi4.g3<-tapply(data$Pi_g3,data$rank,sum)/tapply(data$L4.correct,data$rank,sum)*3
	sd.pi4<-tapply(data$Pi_g3*3/data$L4,data$rank,sd)
	
	TD<-get_TD_bins(data,num.bins, numChr,freq0,g2,site='all')
	TD0<-get_TD_bins(data,num.bins, numChr,freq0,g2,site='0')
	TD4<-get_TD_bins(data,num.bins, numChr,freq0,g2,site='4')
	FWH_norm<-get_FWH_norm_bins(data,num.bins, numChr,freq0,g2, site='all')
	FWH_norm0<-get_FWH_norm_bins(data,num.bins, numChr,freq0,g2, site='0')	
	FWH_norm4<-get_FWH_norm_bins(data,num.bins, numChr,freq0,g2, site='4')	
	
	
	# TD<-tapply(data$D, data$rank, mean, na.rm=T)
	# FWH<-tapply(data$FWH,data$rank,mean,na.rm=T)
	if(is.null(rec.file)){
		res<-data.frame(pi4=pi4.g3,pi0.pi4,pi0,sd.pi4, TD, TD0,TD4, FWH_norm=FWH_norm[,2], FWH_norm0=FWH_norm0[,2], FWH_norm4=FWH_norm4[,2], E=FWH_norm[,1], E0=FWH_norm0[,1],E4=FWH_norm4[,1],wt=wt)
	}else{
		rec<-tapply(data$rec,data$rank,mean,na.rm=T)
		res<-data.frame(pi4=pi4.g3,pi0.pi4,pi0,sd.pi4, TD, TD0,TD4, FWH_norm=FWH_norm[,2], FWH_norm0=FWH_norm0[,2], FWH_norm4=FWH_norm4[,2],E=FWH_norm[,1], E0=FWH_norm0[,1],E4=FWH_norm4[,1],wt=wt,rec)
	}
	if(!is.null(GC)){
		gc3<-tapply(data$GC3,data$rank,mean,na.rm=T)
		res$GC3<-gc3
	}
	if(dfe){
		dir.create(paste('bin',num.bins,sep=''))
		for(j in 1:num.bins){
			geneList.bin<-data$Transcript[data$rank==j]
			total.L0.bin<-round(sum(data$L0.correct[data$rank==j])/3)
			total.L4.bin<-round(sum(data$L4.correct[data$rank==j])/3)
			sfs0.bin<-round(summary(factor(freq0$Freq[freq0$Transcript%in%geneList.bin]))/3)
			sfs4.bin<-summary(factor(freq4$Freq[freq4$Transcript%in%geneList.bin]))
			if(length(sfs0.bin)!= (numChr-1) | length(sfs4.bin)!=(numChr-1)){warnings(paste("Not enough SNP in bin", j))}
				
			if(is.null(base)){
				comment=paste('data',nrow(data),'genes',sep='_')
				dfe.file<-paste(outdir,paste('bin',num.bins,sep=''), paste('dfe_input_bin',j,'.txt',sep=''),sep='/')
			}else{
				comment=paste(base, nrow(data),'genes',sep='_')
				dfe.file<-paste(outdir,'/',paste('bin',num.bins,sep=''),'/', base,'_dfe_input_bin',j,'.txt',sep='')
			}
			cat(paste(1,1,numChr,sep="\t"),"\n",sep='',file=dfe.file)
			cat(sfs4.bin,total.L4.bin,sep="\t",append=T, file=dfe.file)
			cat("\n", sep='',append=T,file=dfe.file)
			cat(sfs0.bin,total.L0.bin,sep="\t",append=T, file=dfe.file)
			cat("\n", sep='',append=T,file=dfe.file)
		}
	}
	return(res)
}






get_TD_bins=function(data, num.bins, numChr, freq0, g2,  site='all'){ ## g2 freq4 subset2
	TD=numeric(num.bins)
	for(i in 1:num.bins){
		#theta_pi4=data$Pi4[data$rank==i]
		theta_pi4=data$Pi_g2[data$rank==i]
		theta_pi0=data$Pi0[data$rank==i]
		#S4=data$SNP4[data$rank==i]
		
		# S0=sum(data$SNP0[data$rank==i])
		transcript_id=as.character(data$Transcript[data$rank==i])
		S0=nrow(freq0[freq0$Transcript%in%transcript_id,])
		S4=nrow(g2[g2$Transcript%in%transcript_id,])
		
		if(site == 'all'){
			freq0_bins=freq0$Freq[freq0$Transcript%in%transcript_id]
			freq4_bins=g2$Freq[g2$Transcript%in%transcript_id]
			TD[i]=add_tajimasD(c(freq4_bins,freq0_bins), S4+S0, numChr)
		}else if(site == '4'){
			freq4_bins=g2$Freq[g2$Transcript%in%transcript_id]
			TD[i]=add_tajimasD(freq4_bins, S4, numChr)
		}else if(site == '0'){
			freq0_bins=freq0$Freq[freq0$Transcript%in%transcript_id]
			TD[i]=add_tajimasD(freq0_bins, S0, numChr)
		}else{
			stop("Wrong code for site value")
		}
	}
	return(TD)
}

add_tajimasD_per_gene=function(theta_pi, S, n){
	a1=0
	a2=0
	for(i in 1:(n-1)){
		a1=a1+1/i
		a2=a2+1/i^2
	}
	
	b1=(n+1)/(n-1)/3
	b2=2*(n^2+n+3)/(9*n*(n-1))
	c1=b1-1/a1
	c2=b2-(n+2)/(a1*n) + a2/a1^2
	e1=c1/a1
	e2=c2/(a1^2+a2)
	D=(theta_pi-S/a1)/sqrt(e1*S+e2*S*(S-1))
	return(D)
}

add_tajimasD=function(freq,S,n){
	num.alleles=freq*n
	theta_pi=sum(num.alleles*(n-num.alleles))/choose(n,2)
	
	a1=0
	a2=0
	for(i in 1:(n-1)){
		a1=a1+1/i
		a2=a2+1/i^2
	}
	b1=(n+1)/(n-1)/3
	b2=2*(n^2+n+3)/(9*n*(n-1))
	c1=b1-1/a1
	c2=b2-(n+2)/(a1*n) + a2/a1^2
	e1=c1/a1
	e2=c2/(a1^2+a2)
	
	D=(theta_pi-S/a1)/sqrt(e1*S+e2*S*(S-1))
	return(D)
}	
	
get_FWH_norm_bins=function(data, num.bins, numChr, freq0,g2, site='all'){
	FWH=matrix(numeric(num.bins*2),ncol=2)
	for(i in 1:num.bins){
# 		theta_pi4=data$Pi4[data$rank==i]
		theta_pi4=data$Pi_g2[data$rank==i]

		theta_pi0=data$Pi0[data$rank==i]
# 		S4=data$SNP4[data$rank==i]
		# S0=sum(data$SNP0[data$rank==i])
		transcript_id=as.character(data$Transcript[data$rank==i])
		S0=nrow(freq0[freq0$Transcript%in%transcript_id,])
		S4=nrow(g2[g2$Transcript%in%transcript_id,])

		if(site=='all'){
			freq0_bins=freq0$Freq[freq0$Transcript%in%transcript_id]
			freq4_bins=g2$Freq[g2$Transcript%in%transcript_id]
			FWH[i,]<-add_fayWu_H_norm(c(freq4_bins,freq0_bins), S4+S0, numChr)
		}else if(site=='4'){
			freq4_bins=g2$Freq[g2$Transcript%in%transcript_id]
			FWH[i,]<-add_fayWu_H_norm(freq4_bins, S4, numChr)
		}else if(site== '0'){
			freq0_bins=freq0$Freq[freq0$Transcript%in%transcript_id]
			FWH[i,]<-add_fayWu_H_norm(freq0_bins, S0, numChr)
		}else{
			stop("Wrong code for site value")
		}
	}
	return(FWH)
}

add_fayWu_H_norm=function(freq, S,n){
	num.alleles=freq*n
	theta_pi=sum(num.alleles*(n-num.alleles))/choose(n,2)
	theta_H=sum(num.alleles^2)/choose(n,2)	
	FWH=(theta_pi-theta_H)
	FWH_norm=(theta_pi-theta_H)/2
	
	## variance
	a1=0
	a2=0
	for(i in 1:(n-1)){
		a1=a1+1/i
		a2=a2+1/i^2
	}
	bn1=a2+1/n^2
	thetaW=S/a1
	thetaW2=S*(S-1)/(a1^2+a2)
	var_H=(n-2)*thetaW/(6*(n-1)) + (18*n^2*(3*n+2)*bn1 - (88*n^3 + 9*n^2-13*n+6))*thetaW2/(9*n*(n-1)^2)
	FWH_norm=FWH_norm/sqrt(var_H)
	
	theta_L=(theta_H+thetaW)/2
	E=theta_L-thetaW
	var_E=(0.5*n/(n-1)-1/a1)*thetaW+(a2/a1^2 + 2*(n/(n-1))^2*a2-2*(n*a2-n+1)/(a1*(n-1))-(3*n+1)/(n-1))*thetaW2
	E=E/sqrt(var_E)
	return(c(E,FWH_norm))
}
	





get_num_pi4_site_hyper=function(num.pi4, L4, num.group=3){ # number of pi4, total length of L4
	num.noPoly=L4-num.pi4
	num.ps1=rhyper(1, num.pi4, num.noPoly, L4/num.group)
	
	num.ps.left=num.pi4-num.ps1
	num.noPoly.left=num.ps.left*(num.noPoly/num.pi4)
	num.ps2=rhyper(1, num.ps.left, num.noPoly.left, L4/num.group)
	
	num.ps3=num.ps.left-num.ps2
	
	res=c(num.ps1,num.ps2,num.ps3)
	names(res)<-c('P1','P2','P3')
	return(res)
}
get_pi_from_freq=function(freq,n){
	minor_alleles=freq*n
	major_alleles=n-minor_alleles
	pi=(minor_alleles*major_alleles)/choose(n,2)
	return(pi)
}

add_pi_by_group=function(data,g, n,index){ ## data == new.data, g==group.freq, n=samples size
	g$Pi4<-(n-g$Freq*n)*g$Freq*n/choose(n,2)
	pi4.genes<-tapply(g$Pi4,factor(g$Transcript),sum)
	# data[,paste('Pi_g',index,sep='')]<-0
	
	id<-names(pi4.genes)
	x<-data.frame(id,pi4.genes)
	var=paste('Pi_g',index,sep='')
	names(x)<-c('Transcript',var)
	
	data1<-data[data$Transcript%in%id,]
	data2<-data[!data$Transcript%in%id,]
	
	data1<-merge(data1,x,by='Transcript')
	data2[,var]<-0
	return(rbind(data1,data2)) ## this step actually create orders for the data by pi !!!! 
}

add_thetaW_by_group=function(data,g, n,index){ ## data == new.data, g==group.freq, n=samples size
	numS<-summary(g$Transcript,maxsum=1e8)
	id<-names(numS)
	a=0
	for(i in 1:(n-1)){
		a=a+1/i
	}
	thetaW=numS/a
	
	x<-data.frame(id,thetaW)
	var=paste('Pi_g',index,sep='')
	names(x)<-c('Transcript',var)
	
	data1<-data[data$Transcript%in%id,]
	data2<-data[!data$Transcript%in%id,]
	
	data1<-merge(data1,x,by='Transcript')
	if(nrow(data2)>0){data2[,var]<-0}
	return(rbind(data1,data2)) ## this step actually create orders for the data by pi !!!! 
}







