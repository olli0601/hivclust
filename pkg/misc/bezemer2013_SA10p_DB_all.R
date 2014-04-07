# The function below is added from "https://github.com/olli0601/abc.n/blob/master/pkg/misc/nabc.startme.R".
my.fade.col<-function(col,alpha=0.5)
{
	return(rgb(col2rgb(col)[1]/255,col2rgb(col)[2]/255,col2rgb(col)[3]/255,alpha))
}
hivc.db.Date2numeric<- function( x )
{
	x	<- as.POSIXlt(x)
	tmp	<- x$year + 1900
	x	<- tmp + round( x$yday / ifelse((tmp%%4==0 & tmp%%100!=0) | tmp%%400==0,366,365), d=3 )
	x	
}


#project.bezemer2013a.figs.v130821<- function()
project.bezemer2013a.figs.v131023_DB<- function()
{
	options(warn=1)
	require(data.table)
	require(RColorBrewer)	 	
	if(0)
	{
		dir.name			<-  "C:/Users/dobezemer/WORK/sequences_20111028/seqB/PAPER/Fig1"
		#signat				<- "130827"
		f.name				<- paste("Fig1",".csv",sep='') #signat
		df					<- read.csv(paste(dir.name,f.name,sep='/'), stringsAsFactors=FALSE)
		
	}
	if(1)
	{
		dir.name			<- "/Users/Oliver/duke/2013_HIV_NL/Bezemer2013_clusters"	
		signat				<- "140402"
		f.name				<- paste("Bezemer2014_clusters_Fig1_",signat,".csv",sep='')
		df					<- read.csv(paste(dir.name,f.name,sep='/'), stringsAsFactors=FALSE)
		f.name				<- paste("Bezemer2014_clusters_R_",signat,".csv",sep='')
		tmp					<- read.csv(paste(dir.name,f.name,sep='/'), stringsAsFactors=FALSE)
		colnames(tmp)[colnames(tmp)=="ClusterName"]<- 'clustername'
		df					<- merge(df, tmp, by='clustername', all.x=1)
		f.name				<- paste("Bezemer2014_clusters_MRCAs91MSM_",signat,".csv",sep='')
		tmp					<- read.csv(paste(dir.name,f.name,sep='/'), stringsAsFactors=FALSE)
		tmp					<- tmp[,c('clustername','MRCA','colour.cross')]
		df					<- merge(df, tmp, by='clustername', all.x=1)
	}
	#	
	CLUSTER.N.THRESHOLD	<- 10
	ylab				<- ""
	yheight				<- 1.5
	cex.points			<- 0.75
	xlab				<- "first date of HIV-1 diagnosis"
	xlab.R				<- "R current"
	outfile				<- "BezemerFig1_tryrecentageDB_all"
	outsignat			<- "140402_or"
	verbose				<- 1
	xlim.space			<- 1.5
	xlim.stretch		<- 2
	
	#deal with missing data and get into reasonable shape
	df[df[,"fpos1"]=="","fpos1"]						<- NA
	df[,"fpos1"]										<- hivc.db.Date2numeric( as.Date(df[,"fpos1"], format="%d/%m/%Y") )
	df[,"MRCA"]											<- gsub("\\s","", df[,"MRCA"])
	df[,"MRCA"]											<- gsub(",",".", df[,"MRCA"],fixed=TRUE)
	#df[df[,"MRCA"]=="","MRCA"]							<- NA
	df[,"MRCA"]											<- as.numeric(df[,"MRCA"])
	df[which(is.na(df[,"recent"])),"recent"]			<- 2
	df[,"recent"]										<- factor(df[,"recent"],levels=c(0,1),labels=c("Unknown","Recent"))
	df[,"patient"]										<- factor(df[,"patient"])
	
	tmp													<- !is.na(df[,"RegionOrigin"]) & (df[,"RegionOrigin"]=="" | df[,"RegionOrigin"]=="U") 		
	df[tmp,"RegionOrigin"]								<- NA
	df[,"RegionOrigin"]									<- factor(df[,"RegionOrigin"])	
	tmp													<- unclass(df[,"RegionOrigin"]	)

	df[which(df[,"diaggroup"]==""),"diaggroup"]			<- NA
	df[,"diaggroup"]									<- factor(df[,"diaggroup"])

	df[which(df[,"agegroup"]==""),"agegroup"]			<- NA
	df[,"agegroup"]										<- factor(df[,"agegroup"])

	df[which(df[,"colour.cross"]=="blue cross"),"colour.cross"]	<- 'blue'
	df[which(df[,"colour.cross"]=="grey cross"),"colour.cross"]	<- 'grey50'
	df[which(df[,"colour.cross"]=="orange"),"colour.cross"]		<- 'orange'
	df[which(df[,"colour.cross"]=="orangeblue"),"colour.cross"]	<- 'pink'

	code												<- seq_along(attr(tmp,"levels"))
	names(code)											<- attr(tmp,"levels")
	print(code)
	tmp[ tmp==code[c("NL")] ]							<- 20
	tmp[ tmp==code[c("NLseAnt")] ]						<- 21
	tmp[ tmp==code[c("Lat")] ]							<- 24
	tmp[ tmp%in%code[c("EUC","EUO","EUW")] ]			<- 23
#	tmp[ tmp%in%code[c("AUS","NAM","OAP","SSA","ZAz")] ]<- 24
	tmp[ tmp%in%code[c("AUS","NAM","OAP","SSA","ZAz","NA")] ]<- 25
#
	tmp[ tmp==code[c("SR")] ]	<- 22
	df[,"RegionOrigin"]									<- factor(tmp, levels=as.character(20:25), labels=c("Netherlands","Dutch Antilles","Suriname","Europe","Latin","Other"))

	

	df[which(df[,"Country"]==""),"Country"]				<- NA
	df[which(df[,"clustername"]==0),"clustername"]		<- NA
	df[,"Country"]										<- factor(df[,"Country"])
	df[which(df[,"hregion"]==""),"hregion"]				<- "unknown"
	df[,"hregion"]										<- factor(df[,"hregion"], levels=c("CU","Ea","NH","No","So","UT","ZH","unknown"))	
	#df[,"HospitalRegion"]								<- factor(df[,"HospitalRegion"], levels=c(1,2,3,4,5,6), labels=c("Amst","N","E","S","W","Curu"))	
	df[which(df[,"sexy"]==""),"sexy"]					<- NA
	df[,"sexy"]											<- factor(df[,"sexy"])
	df[which(df[,"rout"]==""),"rout"]					<- "unknown"
	df[which(df[,"rout"]=="U"),"rout"]					<- "unknown"
	df[which(df[,"rout"]=="CH"),"rout"]					<- "other"
	df[which(df[,"rout"]=="BL"),"rout"]					<- "other"	
	df[,"rout"]											<- factor(df[,"rout"])


	#print(levels(df[,"Country"]))
	#print(levels(df[,"hregion"]))
	#print(levels(df[,"sexy"]))
	#print(levels(df[,"rout"]))
	#print(summary(df[,"fpos1"]))
	#print(levels(df[,"diaggroup"]))
	#print(levels(df[,"agegroup"]))
	df	<- data.table(df, key="clustername")
	setnames(df, c("clustername","hregion"),c("cluster","RegionHospital"))
	str(df)
	#
	#	determine median fpos1 time for cluster and add PLUS determine cex for each cluster and add
	#
	tmp			<- df[,	list(	cluster.PosT=	mean(fpos1, na.rm=T), 
								cluster.cex=	log(nlseq[1]),
								cluster.fNL=	length(which(RegionOrigin=="NL")) / length(RegionOrigin)
								),	by=cluster]
	tmp			<- subset(tmp, !is.na(cluster))
	set(tmp, NULL, "cluster.cex", tmp[,cluster.cex] * 0.5/(max(tmp[,cluster.cex]) - min(tmp[,cluster.cex])) )
	set(tmp, NULL, "cluster.cex", tmp[,cluster.cex]  - max(tmp[,cluster.cex]) + 1)	
	df			<- merge(df, tmp, all.x=1, by="cluster")	
	cat(paste("number of sequences in df, n=", nrow(df)))
	#
	#	remove small clusters and singletons
	#
	singletons			<- subset(df,is.na(cluster))
	cat(paste("number of singletones, n=", nrow(singletons)))
#	clusters.small		<- df[, list(select= (NLseq>1 & nlseq<CLUSTER.N.THRESHOLD)[1]), by= cluster]
	clusters.small		<- df[, list(select= (nlseq>1 & nlseq<CLUSTER.N.THRESHOLD)[1]), by= cluster]	
	clusters.small		<- merge(subset(clusters.small, select), df, by="cluster")
	cat(paste("number seq in small clusters, n=", nrow(clusters.small)))
	df					<- subset(df, !is.na(cluster) & nlseq>=CLUSTER.N.THRESHOLD)
	cat(paste("number seq in large clusters, n=", nrow(df)))
	#
	#	determine if network concentrated in region origin, and where
	#
	clu.reg				<- table( subset(df,RegionHospital!="unknown")[,cluster,RegionHospital] )
	clu.reg				<- clu.reg / matrix(apply(clu.reg, 2, sum), nrow=nrow(clu.reg), ncol=ncol(clu.reg), byrow=1)				
	#ALTERNATIVE DEFINITION
	#clu.regconc	<- apply( clu.reg, 2, function(x)	any(x>0.4) && length(which(x!=0))>2 ) |
	#					apply( clu.reg, 2, function(x)	any(x>0.5) && length(which(x!=0))<=2 )
	clu.regconc			<- apply( clu.reg, 2, function(x)	any(x>0.5) ) 			
	if(verbose)	cat(paste("number of clusters, n=",length(clu.regconc)))
	if(verbose)	cat(paste("number of spatially concentrated clusters, n=",length(which(clu.regconc))))
	tmp					<- rownames(clu.reg)
	clu.regconc2		<- sapply( seq_len(ncol(clu.reg)), function(j)
								{
									ifelse(!clu.regconc[j], NA, tmp[ which.max(clu.reg[,j]) ] )						
								})
	df.clureg			<- data.table(cluster= as.numeric(colnames( clu.reg )), IsRegionConcentrated= clu.regconc, RegionConcentrated=clu.regconc2, key="cluster" )
	df					<- merge(df, df.clureg, all.x=1, by="cluster" )
	#
	#	determine if network concentrated risk group, and which
	#
	clu.reg				<- table( subset(df,rout!="unknown")[,cluster,rout] )
	clu.reg				<- clu.reg / matrix(apply(clu.reg, 2, sum), nrow=nrow(clu.reg), ncol=ncol(clu.reg), byrow=1)				
	clu.regconc			<- apply( clu.reg, 2, function(x)	any(x>0.5) ) 			
	if(verbose)	cat(paste("number of clusters, n=",length(clu.regconc)))
	if(verbose)	cat(paste("number of clusters concentrated by risk group, n=",length(which(clu.regconc))))
	tmp					<- rownames(clu.reg)
	clu.regconc2		<- sapply( seq_len(ncol(clu.reg)), function(j)
			{
				ifelse(!clu.regconc[j], NA, tmp[ which.max(clu.reg[,j]) ] )						
			})
	df.clureg			<- data.table(cluster= as.numeric(colnames( clu.reg )), IsExpGrConcentrated= clu.regconc, ExpGrConcentrated=clu.regconc2, key="cluster" )
	df					<- merge(df, df.clureg, all.x=1, by="cluster" )
	#
	#
	#
	df[,time:=fpos1 ]		
	set(df, which(df[,!IsExpGrConcentrated | is.na(IsExpGrConcentrated)]), "ExpGrConcentrated", "mixed")
	tmp					<- numeric(nrow(df))	
#	tmp2				<- c("DU","HT","MSM","mixed")
	tmp2				<- c("mixed","DU","HT","MSM")
	for( i in seq_along(tmp2))
		tmp[which(df[,ExpGrConcentrated==tmp2[i]])]	<- i
	df[,sort1:=factor(tmp, levels=1:4, labels=tmp2) ]
	df[,sort2:=max(duration,na.rm=1)-duration]
	#	
	xlim				<- range(c(range( df[,fpos1], na.rm=1 )	, range( df[,MRCA], na.rm=1 )))		
	#extract clusters one by one in desired order, with all the information for plotting
	tmp					<- df[,list(sort1=sort1[1], sort2=sort2[1]),by=cluster]
	clusters.sortby1	<- as.numeric( tmp[,sort1] )
	clusters.sortby2	<- as.numeric( tmp[,sort2] )
	clusters.sortby		<- order(clusters.sortby1,clusters.sortby2, decreasing=F)
	clusters			<- tmp[clusters.sortby,cluster]														
	clusters			<- lapply(clusters,function(x)	subset(df,cluster==x,c(time, cluster,rout, RegionOrigin, diaggroup, agegroup, RegionHospital, recent, cluster.cex, sort1, MRCA, colour.cross, MeanR, R025, R975))		)
	ylim				<- c(-50,length(clusters))
	xlim[1]				<- xlim[1]-700/365		
	df.R				<- data.table(idx=seq_along(clusters), cluster=sapply(clusters, function(x) x[, cluster][1]), R.mean=sapply(clusters, function(x) x[, MeanR][1]), R.025=sapply(clusters, function(x) x[, R025][1]), R.975=sapply(clusters, function(x) x[, R975][1]) )
	xlim.R				<- df.R[, range(c(range(R.mean, na.rm=1),range(R.025, na.rm=1),range(R.975, na.rm=1)))]	
	xlim[3]				<- xlim.stretch*xlim.R[2]+xlim[2]+xlim.space
	#
	#	plot by recent
	#
	ncols				<- length(levels(df[,recent]))
	cols				<- c(brewer.pal(9,"Set1"),"grey50")		
	cols				<- sapply(cols, function(x) my.fade.col(x,0.8))[c(9,1,10)]	
	file				<- paste(dir.name,paste(outfile,"_recent_",gsub('/',':',outsignat),".pdf",sep=''),sep='/')
	cat(paste("\nwrite plot to",file))
	pdf(file,width=7,height=12)
	par(mar=c(4,4,0,0))	
	plot(1,1,type='n',bty='n',xlim=xlim[c(1,3)],ylim=ylim,ylab=ylab, xaxt='n',yaxt='n',xlab='')
	axis(1, seq(round(xlim[1]),round(xlim[2]),by=1), labels = TRUE, cex.axis=0.8 )
	mtext(xlab, side=1, line=3, at=mean(xlim[c(1,2)]))
	dummy	<- sapply(seq_along(clusters),function(i)
			{							
				cluster.cex	<- cex.points  * clusters[[i]][1,cluster.cex]
				cluster.ix	<- order(as.numeric(clusters[[i]][,time]))
				cluster.x	<- as.numeric(clusters[[i]][,time])[cluster.ix]										
				cluster.y	<- rep(i,nrow(clusters[[i]]))
				cluster.z	<- as.numeric(clusters[[i]][,recent])[cluster.ix]
				#print(clusters[[i]][,covariate])
				#if(cluster.z[1]!=4)
				#lines( c(cluster.x[1],xlim[2]), rep(i,2), col=cols[ cluster.z[1] ] )					
				points( cluster.x, cluster.y, col=cols[ cluster.z ], pch=19, cex=cluster.cex )
				if(!is.na(clusters[[i]][,MRCA][1]))
					points( clusters[[i]][,MRCA][1], cluster.y[1], col=clusters[[i]][,colour.cross][1], pch=4, cex=cex.points*0.7 )
			})	
	
	c.lines	<- sapply(clusters,function(x) x[1,sort1])
	c.lines	<- c( which(diff(as.numeric(c.lines))!=0), length(clusters) )
	dummy	<- lapply(c(0,c.lines),function(x)		lines( c(xlim[1],xlim[2]), rep(x+0.5,2), lty=3, lwd=0.75, col="grey50")			)
	mtext(text=levels(df[,sort1]), 2, at=c(0,c.lines[-length(c.lines)]) + diff( c(0,c.lines) )/2, las=1, cex=0.7)	
	legend(xlim[1]-600/365,ylim[2]*0.9,bty='n',pt.bg=cols,pch=21,legend=levels(df[,recent]), col=rep("transparent",ncols))			
	#plot small clusters	
	tmp	<- runif(nrow(clusters.small), min=-23,max=-3)
	points( as.numeric(clusters.small[,fpos1]), tmp, col=cols[ as.numeric(clusters.small[,recent]) ], pch=19, cex= cex.points*0.4 )
	mtext(text="patients\nclusters<10", 2, at=mean(c(-23,-3)), las=1, cex=0.7)
	#plot singletons
	tmp	<- runif(nrow(singletons), min=-47,max=-27)
	points( as.numeric(singletons[,fpos1]), tmp, col=cols[ as.numeric(singletons[,recent]) ], pch=19, cex= cex.points*0.4 )
	mtext(text="patients\nnot\nclustering", 2, at=mean(c(-47,-27)), las=1, cex=0.7)
	#	add R current plot
	tmp	<- subset(df.R, !is.na(R.mean))[, {
				lines( xlim[2]+xlim.space+xlim.stretch*c(R.025,R.975), rep(idx,2), lwd=1, lty=1)
				points(xlim[2]+xlim.space+xlim.stretch*R.mean, idx, pch=19, col='black', cex=0.5)
			} ,by=idx]
	lines(rep(xlim[2]+xlim.space+xlim.stretch*1,2), c(ylim), lty=3, lwd=0.75, col="red")	
	axis(1, xlim[2]+xlim.space+xlim.stretch*seq(xlim.R[1],xlim.R[2],by=0.5), labels=seq(xlim.R[1],xlim.R[2],by=0.5), cex.axis=0.8 )
	mtext(xlab.R, side=1, line=3, at=xlim[2]+xlim.space+xlim.stretch*mean(xlim.R))
	dev.off()
	
	#
	#	plot by region origin
	#
	ncols				<- length(levels(df[,RegionOrigin]))
	cols				<- c(brewer.pal(9,"Set1"),"grey50")		
#	cols				<- sapply(cols, function(x) my.fade.col(x,0.8))[c(2,5,3,7,10)]
	cols				<- sapply(cols, function(x) my.fade.col(x,0.8))[c(2,5,3,7,4,10)]	
	file				<- paste(dir.name,paste(outfile,"_RegionOrigin_",gsub('/',':',outsignat),".pdf",sep=''),sep='/')
	cat(paste("\nwrite plot to",file))
	pdf(file,width=7,height=12)
	par(mar=c(4,4,0,0))	
	plot(1,1,type='n',bty='n',xlim=xlim[c(1,3)],ylim=ylim,ylab=ylab, xaxt='n',yaxt='n',xlab='')
	axis(1, seq(round(xlim[1]),round(xlim[2]),by=1), labels = TRUE, cex.axis=0.8 )
	mtext(xlab, side=1, line=3, at=mean(xlim[c(1,2)]))	
	dummy	<- sapply(seq_along(clusters),function(i)
			{							
				cluster.cex	<- cex.points  * clusters[[i]][1,cluster.cex]
				cluster.ix	<- order(as.numeric(clusters[[i]][,time]))
				cluster.x	<- as.numeric(clusters[[i]][,time])[cluster.ix]										
				cluster.y	<- rep(i,nrow(clusters[[i]]))
				cluster.z	<- as.numeric(clusters[[i]][,RegionOrigin])[cluster.ix]
				#print(clusters[[i]][,covariate])
				#if(cluster.z[1]!=4)
				#lines( c(cluster.x[1],xlim[2]), rep(i,2), col=cols[ cluster.z[1] ] )					
				points( cluster.x, cluster.y, col=cols[ cluster.z ], pch=19, cex=cluster.cex )
				if(!is.na(clusters[[i]][,MRCA][1]))
					points( clusters[[i]][,MRCA][1], cluster.y[1], col=clusters[[i]][,colour.cross][1], pch=4, cex=cex.points*0.7 )
			})	
	
	c.lines	<- sapply(clusters,function(x) x[1,sort1])
	c.lines	<- c( which(diff(as.numeric(c.lines))!=0), length(clusters) )
	dummy	<- lapply(c(0,c.lines),function(x)		lines( c(xlim[1],xlim[2]), rep(x+0.5,2), lty=3, lwd=0.75, col="grey50")			)	
	mtext(text=levels(df[,sort1]), 2, at=c(0,c.lines[-length(c.lines)]) + diff( c(0,c.lines) )/2, las=1, cex=0.7)	
	legend(xlim[1]-600/365,ylim[2]*0.9,bty='n',pt.bg=cols,pch=21,legend=levels(df[,RegionOrigin]), col=rep("transparent",ncols))			
	#plot small clusters	
	tmp	<- runif(nrow(clusters.small), min=-23,max=-3)
	points( as.numeric(clusters.small[,fpos1]), tmp, col=cols[ as.numeric(clusters.small[,RegionOrigin]) ], pch=19, cex= cex.points*0.4 )
	mtext(text="patients\nclusters<10", 2, at=mean(c(-23,-3)), las=1, cex=0.7)
	#plot singletons
	tmp	<- runif(nrow(singletons), min=-47,max=-27)
	points( as.numeric(singletons[,fpos1]), tmp, col=cols[ as.numeric(singletons[,RegionOrigin]) ], pch=19, cex= cex.points*0.4 )
	mtext(text="patients\nnot\nclustering", 2, at=mean(c(-47,-27)), las=1, cex=0.7)
	#	add R current plot
	tmp	<- subset(df.R, !is.na(R.mean))[, {
				lines( xlim[2]+xlim.space+xlim.stretch*c(R.025,R.975), rep(idx,2), lwd=1, lty=1)
				points(xlim[2]+xlim.space+xlim.stretch*R.mean, idx, pch=19, col='black', cex=0.5)
			} ,by=idx]
	lines(rep(xlim[2]+xlim.space+xlim.stretch*1,2), c(ylim), lty=3, lwd=0.75, col="red")	
	axis(1, xlim[2]+xlim.space+xlim.stretch*seq(xlim.R[1],xlim.R[2],by=0.5), labels=seq(xlim.R[1],xlim.R[2],by=0.5), cex.axis=0.8 )	
	mtext(xlab.R, side=1, line=3, at=xlim[2]+xlim.space+xlim.stretch*mean(xlim.R))
	dev.off()
	

##

#	plot by diaggroup
	#
	ncols				<- length(levels(df[,diaggroup]))
	cols				<- c(brewer.pal(9,"Set1"),"grey50")		
#	cols				<- sapply(cols, function(x) my.fade.col(x,0.8))[c(2,5,3,7,10, 11, 12)]
	cols				<- sapply(cols, function(x) my.fade.col(x,0.8))[c(2,5,3,7,6,8, 9, 4)]	
	file				<- paste(dir.name,paste(outfile,"_diaggroup_",gsub('/',':',outsignat),".pdf",sep=''),sep='/')
	cat(paste("\nwrite plot to",file))
	pdf(file,width=7,height=12)
	par(mar=c(4,4,0,0))	
	plot(1,1,type='n',bty='n',xlim=xlim[c(1,3)],ylim=ylim,ylab=ylab, xaxt='n',yaxt='n',xlab='')
	axis(1, seq(round(xlim[1]),round(xlim[2]),by=1), labels = TRUE, cex.axis=0.8 )
	mtext(xlab, side=1, line=3, at=mean(xlim[c(1,2)]))
	dummy	<- sapply(seq_along(clusters),function(i)
			{							
				cluster.cex	<- cex.points  * clusters[[i]][1,cluster.cex]
				cluster.ix	<- order(as.numeric(clusters[[i]][,time]))
				cluster.x	<- as.numeric(clusters[[i]][,time])[cluster.ix]										
				cluster.y	<- rep(i,nrow(clusters[[i]]))
				cluster.z	<- as.numeric(clusters[[i]][,diaggroup])[cluster.ix]
				#print(clusters[[i]][,covariate])
				#if(cluster.z[1]!=4)
				#lines( c(cluster.x[1],xlim[2]), rep(i,2), col=cols[ cluster.z[1] ] )					
				points( cluster.x, cluster.y, col=cols[ cluster.z ], pch=19, cex=cluster.cex )
				if(!is.na(clusters[[i]][,MRCA][1]))
					points( clusters[[i]][,MRCA][1], cluster.y[1], col=clusters[[i]][,colour.cross][1], pch=4, cex=cex.points*0.7 )
			})	
	
	c.lines	<- sapply(clusters,function(x) x[1,sort1])
	c.lines	<- c( which(diff(as.numeric(c.lines))!=0), length(clusters) )
	dummy	<- lapply(c(0,c.lines),function(x)		lines( c(xlim[1],xlim[2]), rep(x+0.5,2), lty=3, lwd=0.75, col="grey50")			)	
	mtext(text=levels(df[,sort1]), 2, at=c(0,c.lines[-length(c.lines)]) + diff( c(0,c.lines) )/2, las=1, cex=0.7)	
	legend(xlim[1]-600/365,ylim[2]*0.9,bty='n',pt.bg=cols,pch=21,legend=levels(df[,diaggroup]), col=rep("transparent",ncols))			
	#plot small clusters	
	tmp	<- runif(nrow(clusters.small), min=-23,max=-3)
	points( as.numeric(clusters.small[,fpos1]), tmp, col=cols[ as.numeric(clusters.small[,diaggroup]) ], pch=19, cex= cex.points*0.4 )
	mtext(text="patients\nclusters<10", 2, at=mean(c(-23,-3)), las=1, cex=0.7)
	#plot singletons
	tmp	<- runif(nrow(singletons), min=-47,max=-27)
	points( as.numeric(singletons[,fpos1]), tmp, col=cols[ as.numeric(singletons[,diaggroup]) ], pch=19, cex= cex.points*0.4 )
	mtext(text="patients\nnot\nclustering", 2, at=mean(c(-47,-27)), las=1, cex=0.7)
	#	add R current plot
	tmp	<- subset(df.R, !is.na(R.mean))[, {
				lines( xlim[2]+xlim.space+xlim.stretch*c(R.025,R.975), rep(idx,2), lwd=1, lty=1)
				points(xlim[2]+xlim.space+xlim.stretch*R.mean, idx, pch=19, col='black', cex=0.5)
			} ,by=idx]
	lines(rep(xlim[2]+xlim.space+xlim.stretch*1,2), c(ylim), lty=3, lwd=0.75, col="red")	
	axis(1, xlim[2]+xlim.space+xlim.stretch*seq(xlim.R[1],xlim.R[2],by=0.5), labels=seq(xlim.R[1],xlim.R[2],by=0.5), cex.axis=0.8 )	
	mtext(xlab.R, side=1, line=3, at=xlim[2]+xlim.space+xlim.stretch*mean(xlim.R))	
	dev.off()


#	plot by agegroup
	#
	ncols				<- length(levels(df[,agegroup]))
	cols				<- c(brewer.pal(9,"Set1"),"grey50")		
#	cols				<- sapply(cols, function(x) my.fade.col(x,0.8))[c(2,5,3,7,10, 11, 12)]
	cols				<- sapply(cols, function(x) my.fade.col(x,0.8))[c(2,5,3,7,6,8, 9, 4)]	
	file				<- paste(dir.name,paste(outfile,"_agegroup_",gsub('/',':',outsignat),".pdf",sep=''),sep='/')
	cat(paste("\nwrite plot to",file))
	pdf(file,width=7,height=12)
	par(mar=c(4,4,0,0))	
	plot(1,1,type='n',bty='n',xlim=xlim[c(1,3)],ylim=ylim,ylab=ylab, xaxt='n',yaxt='n',xlab='')
	axis(1, seq(round(xlim[1]),round(xlim[2]),by=1), labels = TRUE, cex.axis=0.8 )
	mtext(xlab, side=1, line=3, at=mean(xlim[c(1,2)]))
	dummy	<- sapply(seq_along(clusters),function(i)
			{							
				cluster.cex	<- cex.points  * clusters[[i]][1,cluster.cex]
				cluster.ix	<- order(as.numeric(clusters[[i]][,time]))
				cluster.x	<- as.numeric(clusters[[i]][,time])[cluster.ix]										
				cluster.y	<- rep(i,nrow(clusters[[i]]))
				cluster.z	<- as.numeric(clusters[[i]][,agegroup])[cluster.ix]
				#print(clusters[[i]][,covariate])
				#if(cluster.z[1]!=4)
				#lines( c(cluster.x[1],xlim[2]), rep(i,2), col=cols[ cluster.z[1] ] )					
				points( cluster.x, cluster.y, col=cols[ cluster.z ], pch=19, cex=cluster.cex )
				if(!is.na(clusters[[i]][,MRCA][1]))
					points( clusters[[i]][,MRCA][1], cluster.y[1], col=clusters[[i]][,colour.cross][1], pch=4, cex=cex.points*0.7 )
			})	
	
	c.lines	<- sapply(clusters,function(x) x[1,sort1])
	c.lines	<- c( which(diff(as.numeric(c.lines))!=0), length(clusters) )
	dummy	<- lapply(c(0,c.lines),function(x)		lines( c(xlim[1],xlim[2]), rep(x+0.5,2), lty=3, lwd=0.75, col="grey50")			)	
	mtext(text=levels(df[,sort1]), 2, at=c(0,c.lines[-length(c.lines)]) + diff( c(0,c.lines) )/2, las=1, cex=0.7)	
	legend(xlim[1]-600/365,ylim[2]*0.9,bty='n',pt.bg=cols,pch=21,legend=levels(df[,agegroup]), col=rep("transparent",ncols))			
	#plot small clusters	
	tmp	<- runif(nrow(clusters.small), min=-23,max=-3)
	points( as.numeric(clusters.small[,fpos1]), tmp, col=cols[ as.numeric(clusters.small[,agegroup]) ], pch=19, cex= cex.points*0.4 )
	mtext(text="patients\nclusters<10", 2, at=mean(c(-23,-3)), las=1, cex=0.7)
	#plot singletons
	tmp	<- runif(nrow(singletons), min=-47,max=-27)
	points( as.numeric(singletons[,fpos1]), tmp, col=cols[ as.numeric(singletons[,agegroup]) ], pch=19, cex= cex.points*0.4 )
	mtext(text="patients\nnot\nclustering", 2, at=mean(c(-47,-27)), las=1, cex=0.7)	
	#	add R current plot
	tmp	<- subset(df.R, !is.na(R.mean))[, {
				lines( xlim[2]+xlim.space+xlim.stretch*c(R.025,R.975), rep(idx,2), lwd=1, lty=1)
				points(xlim[2]+xlim.space+xlim.stretch*R.mean, idx, pch=19, col='black', cex=0.5)
			} ,by=idx]
	lines(rep(xlim[2]+xlim.space+xlim.stretch*1,2), c(ylim), lty=3, lwd=0.75, col="red")	
	axis(1, xlim[2]+xlim.space+xlim.stretch*seq(xlim.R[1],xlim.R[2],by=0.5), labels=seq(xlim.R[1],xlim.R[2],by=0.5), cex.axis=0.8 )	
	mtext(xlab.R, side=1, line=3, at=xlim[2]+xlim.space+xlim.stretch*mean(xlim.R))		
	dev.off()


	#
	#	plot by exposure group
	#	
	ncols				<- length(levels(df[,rout]))
	cols				<- c(brewer.pal(ncols-1,"Set1"),"grey50")		
	cols				<- sapply(cols, function(x) my.fade.col(x,0.7))[c(1,3,2,4,5)]	
	file				<- paste(dir.name,paste(outfile,"_rout_",gsub('/',':',outsignat),".pdf",sep=''),sep='/')
	cat(paste("\nwrite plot to",file))
	pdf(file,width=7,height=12)
	par(mar=c(4,4,0,0))	
	plot(1,1,type='n',bty='n',xlim=xlim[c(1,3)],ylim=ylim,ylab=ylab, xaxt='n',yaxt='n',xlab='')
	axis(1, seq(round(xlim[1]),round(xlim[2]),by=1), labels = TRUE, cex.axis=0.8 )
	mtext(xlab, side=1, line=3, at=mean(xlim[c(1,2)]))
	dummy	<- sapply(seq_along(clusters),function(i)
			{							
				cluster.cex	<- cex.points  * clusters[[i]][1,cluster.cex]
				cluster.ix	<- order(as.numeric(clusters[[i]][,time]))
				cluster.x	<- as.numeric(clusters[[i]][,time])[cluster.ix]										
				cluster.y	<- rep(i,nrow(clusters[[i]]))
				cluster.z	<- as.numeric(clusters[[i]][,rout])[cluster.ix]
				#print(clusters[[i]][,covariate])
				#if(cluster.z[1]!=4)
				#lines( c(cluster.x[1],xlim[2]), rep(i,2), col=cols[ cluster.z[1] ] )					
				points( cluster.x, cluster.y, col=cols[ cluster.z ], pch=19, cex=cluster.cex )
				if(!is.na(clusters[[i]][,MRCA][1]))
					points( clusters[[i]][,MRCA][1], cluster.y[1], col=clusters[[i]][,colour.cross][1], pch=4, cex=cex.points*0.7 )
			})	
	c.lines	<- sapply(clusters,function(x) x[1,sort1])
	c.lines	<- c( which(diff(as.numeric(c.lines))!=0), length(clusters) )
	dummy	<- lapply(c(0,c.lines),function(x)		lines( c(xlim[1],xlim[2]), rep(x+0.5,2), lty=3, lwd=0.75, col="grey50")			)	
	mtext(text=levels(df[,sort1]), 2, at=c(0,c.lines[-length(c.lines)]) + diff( c(0,c.lines) )/2, las=1, cex=0.7)	
	legend(xlim[1]-600/365,ylim[2]*0.9,bty='n',pt.bg=cols,pch=21,legend=levels(df[,rout]), col=rep("transparent",ncols))			
	#plot small clusters	
	tmp	<- runif(nrow(clusters.small), min=-23,max=-3)
	points( as.numeric(clusters.small[,fpos1]), tmp, col=cols[ as.numeric(clusters.small[,rout]) ], pch=19, cex= cex.points*0.4 )
	mtext(text="patients\nclusters<10", 2, at=mean(c(-23,-3)), las=1, cex=0.7)
	#plot singletons
	tmp	<- runif(nrow(singletons), min=-47,max=-27)
	points( as.numeric(singletons[,fpos1]), tmp, col=cols[ as.numeric(singletons[,rout]) ], pch=19, cex= cex.points*0.4 )
	mtext(text="patients\nnot\nclustering", 2, at=mean(c(-47,-27)), las=1, cex=0.7)
	#	add R current plot
	tmp	<- subset(df.R, !is.na(R.mean))[, {
				lines( xlim[2]+xlim.space+xlim.stretch*c(R.025,R.975), rep(idx,2), lwd=1, lty=1)
				points(xlim[2]+xlim.space+xlim.stretch*R.mean, idx, pch=19, col='black', cex=0.5)
			} ,by=idx]
	lines(rep(xlim[2]+xlim.space+xlim.stretch*1,2), c(ylim), lty=3, lwd=0.75, col="red")	
	axis(1, xlim[2]+xlim.space+xlim.stretch*seq(xlim.R[1],xlim.R[2],by=0.5), labels=seq(xlim.R[1],xlim.R[2],by=0.5), cex.axis=0.8 )	
	mtext(xlab.R, side=1, line=3, at=xlim[2]+xlim.space+xlim.stretch*mean(xlim.R))			
	dev.off()	
	
	
	
	
	
	#
	#	plot by RegionHospital
	#
	df[,time:=fpos1 ]		
	set(df, which(df[,!IsRegionConcentrated | is.na(IsRegionConcentrated)]), "RegionConcentrated", "mixed")
	tmp					<- numeric(nrow(df))	
	tmp2				<- c("NH","CU","Ea","No","So","UT","ZH","mixed")
	for( i in seq_along(tmp2))
		tmp[which(df[,RegionConcentrated==tmp2[i]])]	<- i
	df[,sort1:=factor(tmp, levels=1:8, labels=tmp2) ]
	df[,sort2:=max(duration,na.rm=1)-duration]
	#	
	xlim				<- range( df[,fpos1], na.rm=1 )		
	#extract clusters one by one in desired order, with all the information for plotting
	tmp					<- df[,list(sort1=sort1[1], sort2=sort2[1]),by=cluster]
	clusters.sortby1	<- as.numeric( tmp[,sort1] )
	clusters.sortby2	<- as.numeric( tmp[,sort2] )
	clusters.sortby		<- order(clusters.sortby1,clusters.sortby2, decreasing=F)
	clusters			<- tmp[clusters.sortby,cluster]		
	clusters			<- lapply(clusters,function(x)	subset(df,cluster==x,c(time, cluster,rout, RegionOrigin, diaggroup, agegroup, RegionHospital, recent, cluster.cex, sort1, MRCA, colour.cross, MeanR, R025, R975))		)
	ylim				<- c(-50,length(clusters))
	xlim[1]				<- xlim[1]-700/365		
	df.R				<- data.table(idx=seq_along(clusters), cluster=sapply(clusters, function(x) x[, cluster][1]), R.mean=sapply(clusters, function(x) x[, MeanR][1]), R.025=sapply(clusters, function(x) x[, R025][1]), R.975=sapply(clusters, function(x) x[, R975][1]) )
	xlim.R				<- df.R[, range(c(range(R.mean, na.rm=1),range(R.025, na.rm=1),range(R.975, na.rm=1)))]	
	xlim[3]				<- xlim.stretch*xlim.R[2]+xlim[2]+xlim.space
	
	ncols				<- length(levels(df[,RegionHospital]))
	cols				<- c(brewer.pal(9,"Set1"),"grey50")		
	cols				<- sapply(cols, function(x) my.fade.col(x,0.8))[c(1,3,2,4,5,6,7,10)]	
	file				<- paste(dir.name,paste(outfile,"_RegionHospital_",gsub('/',':',outsignat),".pdf",sep=''),sep='/')
	cat(paste("\nwrite plot to",file))
	pdf(file,width=7,height=12)
	par(mar=c(4,4,0,0))	
	plot(1,1,type='n',bty='n',xlim=xlim[c(1,3)],ylim=ylim,ylab=ylab, xaxt='n',yaxt='n',xlab='')
	axis(1, seq(round(xlim[1]),round(xlim[2]),by=1), labels = TRUE, cex.axis=0.8 )
	mtext(xlab, side=1, line=3, at=mean(xlim[c(1,2)]))
	dummy	<- sapply(seq_along(clusters),function(i)
			{							
				cluster.cex	<- cex.points  * clusters[[i]][1,cluster.cex]
				cluster.ix	<- order(as.numeric(clusters[[i]][,time]))
				cluster.x	<- as.numeric(clusters[[i]][,time])[cluster.ix]										
				cluster.y	<- rep(i,nrow(clusters[[i]]))
				cluster.z	<- as.numeric(clusters[[i]][,RegionHospital])[cluster.ix]
				#print(clusters[[i]][,covariate])
				#if(cluster.z[1]!=4)
				#lines( c(cluster.x[1],xlim[2]), rep(i,2), col=cols[ cluster.z[1] ] )					
				points( cluster.x, cluster.y, col=cols[ cluster.z ], pch=19, cex=cluster.cex )
				if(!is.na(clusters[[i]][,MRCA][1]))
					points( clusters[[i]][,MRCA][1], cluster.y[1], col=clusters[[i]][,colour.cross][1], pch=4, cex=cex.points*0.7 )
			})	
	
	tmp		<- sapply(clusters,function(x) x[1,sort1])
	c.lines	<- c( which(diff(as.numeric(tmp))!=0), length(clusters) )
	dummy	<- lapply(c(0,c.lines),function(x)		lines( c(xlim[1],xlim[2]), rep(x+0.5,2), lty=3, lwd=0.75, col="grey50")			)
	mtext(text=unique(tmp), 2, at=c(0,c.lines[-length(c.lines)]) + diff( c(0,c.lines) )/2, las=1, cex=0.5)	
	legend(xlim[1]-600/365,ylim[2]*1,bty='n',pt.bg=cols,pch=21,legend=levels(df[,RegionHospital]), col=rep("transparent",ncols))			
	#plot small clusters	
	tmp	<- runif(nrow(clusters.small), min=-23,max=-3)
	points( as.numeric(clusters.small[,fpos1]), tmp, col=cols[ as.numeric(clusters.small[,RegionHospital]) ], pch=19, cex= cex.points*0.4 )
	mtext(text="patients\nclusters<10", 2, at=mean(c(-23,-3)), las=1, cex=0.7)
	#plot singletons
	tmp	<- runif(nrow(singletons), min=-47,max=-27)
	points( as.numeric(singletons[,fpos1]), tmp, col=cols[ as.numeric(singletons[,RegionHospital]) ], pch=19, cex= cex.points*0.4 )
	mtext(text="patients\nnot\nclustering", 2, at=mean(c(-47,-27)), las=1, cex=0.7)
	#	add R current plot
	tmp	<- subset(df.R, !is.na(R.mean))[, {
				lines( xlim[2]+xlim.space+xlim.stretch*c(R.025,R.975), rep(idx,2), lwd=1, lty=1)
				points(xlim[2]+xlim.space+xlim.stretch*R.mean, idx, pch=19, col='black', cex=0.5)
			} ,by=idx]
	lines(rep(xlim[2]+xlim.space+xlim.stretch*1,2), c(ylim), lty=3, lwd=0.75, col="red")	
	axis(1, xlim[2]+xlim.space+xlim.stretch*seq(xlim.R[1],xlim.R[2],by=0.5), labels=seq(xlim.R[1],xlim.R[2],by=0.5), cex.axis=0.8 )	
	mtext(xlab.R, side=1, line=3, at=xlim[2]+xlim.space+xlim.stretch*mean(xlim.R))				
	dev.off()
	
}


#project.bezemer2013b.rates<- function()
#{
#	project.bezemer2013b.rates.v131023()	
#}
#project.bezemer2013a.figs.v131023_DB()