project.bezemer2013a.figs.v130715<- function()
{
	require(data.table)
	require(RColorBrewer)	 
	verbose			<- 1
	dir.name		<- "/Users/Oliver/duke/2013_HIV_NL/Bezemer2013_clusters"	
	signat			<- "130715"
	infile			<- "Bezemer2014_clusters_"
	f.name			<- paste(infile,signat,".csv",sep='')
	df				<- read.csv(paste(dir.name,f.name,sep='/'), stringsAsFactors=FALSE)
	
	#deal with missing data and get into reasonable shape
	#@TODO Daniela, see if I interpret all the fields correctly, especially what is NA
	df[df[,"fpos1"]=="","fpos1"]						<- NA
	df[,"fpos1"]										<- as.Date(df[,"fpos1"], format="%d/%m/%Y")	
	df[,"patient"]										<- factor(df[,"patient"])
	tmp													<- !is.na(df[,"RegionOrigin"]) & (df[,"RegionOrigin"]=="" | df[,"RegionOrigin"]=="U") 		
	df[tmp,"RegionOrigin"]								<- NA
	df[,"RegionOrigin"]									<- factor(df[,"RegionOrigin"])
	
	tmp													<- unclass(df[,"RegionOrigin"]	)
	code												<- seq_along(attr(tmp,"levels"))
	names(code)											<- attr(tmp,"levels")
	print(code)
	tmp[ tmp==code[c("NL")] ]							<- 20
	tmp[ tmp==code[c("NLseAnt")] ]						<- 21
	tmp[ tmp==code[c("Lat")] ]							<- 22
	tmp[ tmp%in%code[c("EUC","EUO","EUW")] ]			<- 23
	tmp[ tmp%in%code[c("AUS","NAM","OAP","SSA","ZAz")] ]<- 24
	df[,"RegionOrigin"]									<- factor(tmp, levels=as.character(20:24), labels=c("NL","NLseAnt","Latin","Europe","other"))
	
	
	df[which(df[,"Country"]==""),"Country"]				<- NA
	df[which(df[,"clustername"]==0),"clustername"]		<- NA
	df[,"Country"]										<- factor(df[,"Country"])
	df[which(df[,"hregion"]==""),"hregion"]				<- "unknown"
	df[,"hregion"]										<- factor(df[,"hregion"], levels=c("CU","Ea","NH","No","So","UT","ZH","unknown"))	
	df[,"HospitalRegion"]								<- factor(df[,"HospitalRegion"], levels=c(1,2,3,4,5,6), labels=c("Amst","N","E","S","W","Curu"))	
	df[which(df[,"sexy"]==""),"sexy"]					<- NA
	df[,"sexy"]											<- factor(df[,"sexy"])
	df[which(df[,"rout"]==""),"rout"]					<- "unknown/other"
	df[which(df[,"rout"]=="CH"),"rout"]					<- "unknown/other"
	df[which(df[,"rout"]=="BL"),"rout"]					<- "unknown/other"
	df[which(df[,"rout"]=="U"),"rout"]					<- "unknown/other"
	df[,"rout"]											<- factor(df[,"rout"])
	#print(levels(df[,"Country"]))
	#print(levels(df[,"hregion"]))
	#print(levels(df[,"sexy"]))
	#print(levels(df[,"rout"]))
	#print(summary(df[,"fpos1"]))
	df	<- data.table(df, key="clustername")
	setnames(df, c("clustername","HospitalRegion"),c("cluster","RegionHospital"))
	str(df)
		
	#determine if network concentrated in region origin, and where
	clu.reg		<- table( df[,cluster,RegionOrigin] )
	clu.reg		<- clu.reg / matrix(apply(clu.reg, 2, sum), nrow=nrow(clu.reg), ncol=ncol(clu.reg), byrow=1)				
	#ALTERNATIVE DEFINITION
	#clu.regconc	<- apply( clu.reg, 2, function(x)	any(x>0.4) && length(which(x!=0))>2 ) |
	#					apply( clu.reg, 2, function(x)	any(x>0.5) && length(which(x!=0))<=2 )
	clu.regconc	<- apply( clu.reg, 2, function(x)	any(x>0.5) ) 			
	if(verbose)	cat(paste("number of clusters, n=",length(clu.regconc)))
	if(verbose)	cat(paste("number of spatially concentrated clusters, n=",length(which(clu.regconc))))
	tmp			<- rownames(clu.reg)
	clu.regconc2<- sapply( seq_len(ncol(clu.reg)), function(j)
			{
				ifelse(!clu.regconc[j], NA, tmp[ which.max(clu.reg[,j]) ] )						
			})
	df.clureg	<- data.table(cluster= as.numeric(colnames( clu.reg )), IsRegionConcentrated= clu.regconc, RegionConcentrated=clu.regconc2, key="cluster" )
	#determine median fpos1 time for cluster and add
	#PLUS determine cex for each cluster and add
	tmp			<- df[,list(cluster.PosT=mean(fpos1, na.rm=T), cluster.cex=log(sequencesperCluster[1])),by=cluster]
	set(tmp, NULL, "cluster.cex", tmp[,cluster.cex] * 0.5/(max(tmp[,cluster.cex]) - min(tmp[,cluster.cex])) )
	set(tmp, NULL, "cluster.cex", tmp[,cluster.cex]  - max(tmp[,cluster.cex]) + 1)	
	df.clureg	<- merge(df.clureg, tmp, all.x=1, by="cluster")	
	#merge all
	df.main		<- merge( df, df.clureg, all.x=1, by="cluster" )
	
	#remove non-clustering sequences
	df.main		<- subset(df.main, !is.na(cluster))	
	
	#print main statistics of networks by region of origin
	print(table(df.main[,RegionConcentrated]))
	
	CLUSTER.N.THRESHOLD	<- 10
	ylab				<- "transmission network ID"
	yheight				<- 1.5
	cex.points			<- 0.5	
	
	#plot by sort1 and sort2	
	df					<- df.main
	df[,time:=fpos1 ]		
	set(df, which(df[,!IsRegionConcentrated]), "RegionConcentrated", "Bridge")
	tmp					<- numeric(nrow(df))	
	tmp2				<- c("NL","Europe","Latin","NLseAnt","other","Bridge")
	for( i in seq_along(tmp2))
		tmp[which(df[,RegionConcentrated==tmp2[i]])]	<- i
	df[,sort1:=factor(tmp, levels=1:6, labels=tmp2) ]
	df[,sort3:=rout ]
	set(df, which(df[,is.na(sort3)]),"sort3","unknown/other")
	df[,sort2:=cluster.PosT]
	
	tmp					<- numeric(nrow(df))		
	tmp2				<- c("NL","Europe","Latin","NLseAnt","other")
	for( i in seq_along(tmp2))
		tmp[which(df[,RegionOrigin==tmp2[i]])]	<- i
	df[,covariate:=factor(tmp, levels=1:5, labels=tmp2) ]
	xlab				<- "first pos date"
	xlim				<- range( df[,fpos1], na.rm=1 )
	ylab				<- ""#paste("clusters BS=",thresh.bs*100," BRL.med=",thresh.brl*100," version=",gsub('/',':',signat),sep="")	
	
	#extract clusters one by one in desired order, with all the information for plotting
	tmp					<- df[,list(sort1=sort1[1], sort2=sort2[1], sort3= as.numeric(length(which(sort3=="MSM"))/length(sort3) > 0.7 )[1]),by=cluster]
	clusters.sortby1	<- as.numeric( tmp[,sort1] )
	clusters.sortby2	<- as.numeric( tmp[,sort2] )
	clusters.sortby3	<- as.numeric( tmp[,sort3] )	
	clusters.sortby		<- order(clusters.sortby1,clusters.sortby3,clusters.sortby2, decreasing=F)
	clusters			<- tmp[clusters.sortby,cluster]														
	clusters			<- lapply(clusters,function(x)	subset(df,cluster==x,c(time,covariate,rout, RegionOrigin, cluster.cex,sort1))		)
	ylim				<- c(-2,length(clusters))
	xlim[1]				<- xlim[1]-700
	
	if(1)	#plot hospital region in colors to justify ordering of plot
	{
		ncols				<- length(levels(df[,covariate]))
		cols				<- c("darkblue","Cyan","firebrick1","darkorange","grey50")	#brewer.pal(ncols,"Set1")#[c(1,6,3,4,5,2)]		
		cols				<- sapply(cols, function(x) my.fade.col(x,0.7))	
		
		
		file				<- paste(dir.name,paste(infile,"RegionOrigin_",gsub('/',':',signat),".pdf",sep=''),sep='/')
		cat(paste("\nwrite plot to",file))
		pdf(file,width=7,height=14)
		plot(1,1,type='n',bty='n',xlim=xlim,ylim=ylim,ylab=ylab, xaxt='n',yaxt='n',xlab=xlab)
		axis.Date(1, seq.Date(xlim[1],xlim[2],by="year"), labels = TRUE )
		dummy	<- sapply(seq_along(clusters),function(i)
				{							
					cluster.cex	<- cex.points  * clusters[[i]][1,cluster.cex]
					cluster.ix	<- order(as.numeric(clusters[[i]][,time]))
					cluster.x	<- as.numeric(clusters[[i]][,time])[cluster.ix]										
					cluster.y	<- rep(i,nrow(clusters[[i]]))
					cluster.z	<- as.numeric(clusters[[i]][,covariate])[cluster.ix]
					#print(clusters[[i]][,covariate])
					#if(cluster.z[1]!=4)
					#lines( c(cluster.x[1],xlim[2]), rep(i,2), col=cols[ cluster.z[1] ] )					
					points( cluster.x, cluster.y, col=cols[ cluster.z ], pch=19, cex=cluster.cex )										
				})	
		c.lines	<- sapply(clusters,function(x) x[1,sort1])
		c.lines	<- which(diff(as.numeric(c.lines))!=0)
		lapply(c.lines,function(x)		abline(h=x+0.5, lty=3, lwd=0.75, col="grey50")			)
		
		mtext(text=levels(df[,covariate]), 2, at=c(0,c.lines[-length(c.lines)]) + diff( c(0,c.lines) )/2, las=1)
		legend("topleft",bty='n',pt.bg=cols,pch=21,legend=levels(df[,covariate]), col=rep("transparent",ncols))				
		dev.off()	
	}
	if(1)	#plot exposure group in colors
	{
		ncols				<- length(levels(df[,rout]))
		cols				<- c(brewer.pal(ncols-1,"Set1"),"grey50")		
		cols				<- sapply(cols, function(x) my.fade.col(x,0.7))[c(1,3,2,4)]	
		
		file				<- paste(dir.name,paste(infile,"rout_",gsub('/',':',signat),".pdf",sep=''),sep='/')
		cat(paste("\nwrite plot to",file))
		pdf(file,width=7,height=14)
		plot(1,1,type='n',bty='n',xlim=xlim,ylim=ylim,ylab=ylab, xaxt='n',yaxt='n',xlab=xlab)
		axis.Date(1, seq.Date(xlim[1],xlim[2],by="year"), labels = TRUE )
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
				})	
		c.lines	<- sapply(clusters,function(x) x[1,sort1])
		c.lines	<- which(diff(as.numeric(c.lines))!=0)
		lapply(c.lines,function(x)		abline(h=x+0.5, lty=3, lwd=0.75, col="grey50")			)
		
		mtext(text=levels(df[,covariate]), 2, at=c(0,c.lines[-length(c.lines)]) + diff( c(0,c.lines) )/2, las=1)
		legend("topleft",bty='n',pt.bg=cols,pch=21,legend=levels(df[,rout]), col=rep("transparent",ncols))				
		dev.off()	
	}
}

project.bezemer2013a.figs.v130714<- function()
{
	require(data.table)
	require(RColorBrewer)	 
	verbose			<- 1
	dir.name		<- "/Users/Oliver/duke/2013_HIV_NL/Bezemer2013_clusters"	
	signat			<- "130714"
	infile			<- "Bezemer2014_clusters_"
	f.name			<- paste(infile,signat,".csv",sep='')
	df				<- read.csv(paste(dir.name,f.name,sep='/'), stringsAsFactors=FALSE)
	
	#deal with missing data and get into reasonable shape
	#@TODO Daniela, see if I interpret all the fields correctly, especially what is NA
	df[df[,"fpos1"]=="","fpos1"]						<- NA
	df[,"fpos1"]										<- as.Date(df[,"fpos1"], format="%d/%m/%Y")	
	df[,"patient"]										<- factor(df[,"patient"])
	tmp													<- !is.na(df[,"RegionOrigin"]) & (df[,"RegionOrigin"]=="" | df[,"RegionOrigin"]=="U") 		
	df[tmp,"RegionOrigin"]								<- NA
	df[,"RegionOrigin"]									<- factor(df[,"RegionOrigin"])
		
	tmp													<- unclass(df[,"RegionOrigin"]	)
	code												<- seq_along(attr(tmp,"levels"))
	names(code)											<- attr(tmp,"levels")
	print(code)
	tmp[ tmp==code[c("NL")] ]							<- 20
	tmp[ tmp==code[c("NLseAnt")] ]						<- 21
	tmp[ tmp==code[c("Lat")] ]							<- 22
	tmp[ tmp%in%code[c("EUC","EUO","EUW")] ]			<- 23
	tmp[ tmp%in%code[c("AUS","NAM","OAP","SSA","ZAz")] ]<- 24
	df[,"RegionOrigin"]									<- factor(tmp, levels=as.character(20:24), labels=c("NL","NLseAnt","Latin","Europe","other"))
	
	
	df[which(df[,"Country"]==""),"Country"]				<- NA
	df[which(df[,"clustername"]==0),"clustername"]		<- NA
	df[,"Country"]										<- factor(df[,"Country"])
	df[which(df[,"hregion"]==""),"hregion"]				<- "unknown"
	df[,"hregion"]										<- factor(df[,"hregion"], levels=c("CU","Ea","NH","No","So","UT","ZH","unknown"))	
	df[,"HospitalRegion"]								<- factor(df[,"HospitalRegion"], levels=c(1,2,3,4,5,6), labels=c("Amst","N","E","S","W","Curu"))	
	df[which(df[,"sexy"]==""),"sexy"]					<- NA
	df[,"sexy"]											<- factor(df[,"sexy"])
	df[which(df[,"rout"]==""),"rout"]					<- "unknown/other"
	df[which(df[,"rout"]=="CH"),"rout"]					<- "unknown/other"
	df[which(df[,"rout"]=="BL"),"rout"]					<- "unknown/other"
	df[which(df[,"rout"]=="U"),"rout"]					<- "unknown/other"
	df[,"rout"]											<- factor(df[,"rout"])
	#print(levels(df[,"Country"]))
	#print(levels(df[,"hregion"]))
	#print(levels(df[,"sexy"]))
	#print(levels(df[,"rout"]))
	#print(summary(df[,"fpos1"]))
	df	<- data.table(df, key="clustername")
	setnames(df, c("clustername","HospitalRegion"),c("cluster","RegionHospital"))
	str(df)
	
	#evaluate distribution by "Country" and collect countries from which there is above level migration to NL
	if(0)	 
	{
		tmp<- sort(table(df[,Country]),decreasing=TRUE)
		barplot(tmp, cex.names=.75)
		tmp<- rev(tmp)/sum(tmp)
		levels<- c(0.01,0.05)
		countries.migration.abovelevel	<- lapply(levels,function(level)
				{
					names(tmp)[cumsum(tmp)>level]
				})
		names(countries.migration.abovelevel)<- levels							
		dput(countries.migration.abovelevel)				
	}	
	#determine if cluster concentrated in region, and where
	clu.reg		<- table( df[,cluster,RegionHospital] )
	clu.reg		<- clu.reg / matrix(apply(clu.reg, 2, sum), nrow=nrow(clu.reg), ncol=ncol(clu.reg), byrow=1)				
	clu.regconc	<- apply( clu.reg, 2, function(x)	any(x>0.4) && length(which(x!=0))>2 ) |
					apply( clu.reg, 2, function(x)	any(x>0.5) && length(which(x!=0))<=2 )
	if(verbose)	cat(paste("number of clusters, n=",length(clu.regconc)))
	if(verbose)	cat(paste("number of spatially concentrated clusters, n=",length(which(clu.regconc))))
	tmp			<- rownames(clu.reg)
	clu.regconc2<- sapply( seq_len(ncol(clu.reg)), function(j)
			{
				ifelse(!clu.regconc[j], NA, tmp[ which.max(clu.reg[,j]) ] )						
			})
	df.clureg	<- data.table(cluster= as.numeric(colnames( clu.reg )), IsRegionConcentrated= clu.regconc, RegionConcentrated=clu.regconc2, key="cluster" )
	#determine median fpos1 time for cluster and add
	#PLUS determine cex for each cluster and add
	tmp			<- df[,list(cluster.PosT=mean(fpos1, na.rm=T), cluster.cex=log(sequencesperCluster[1])),by=cluster]
	set(tmp, NULL, "cluster.cex", tmp[,cluster.cex] * 0.5/(max(tmp[,cluster.cex]) - min(tmp[,cluster.cex])) )
	set(tmp, NULL, "cluster.cex", tmp[,cluster.cex]  - max(tmp[,cluster.cex]) + 1)	
	df.clureg	<- merge(df.clureg, tmp, all.x=1, by="cluster")	
	#merge all
	df.main		<- merge( df, df.clureg, all.x=1, by="cluster" )
	
	#remove non-clustering sequences
	df.main		<- subset(df.main, !is.na(cluster))	

	CLUSTER.N.THRESHOLD	<- 10
	ylab				<- "transmission network ID"
	yheight				<- 1.5
	cex.points			<- 0.5	
	
	#plot by sort1 and sort2	
	df					<- df.main
	df[,time:=fpos1 ]		
	set(df, which(df[,!IsRegionConcentrated]), "RegionConcentrated", "Bridge")
	tmp					<- numeric(nrow(df))
	tmp2				<- c("Curu","N","E","S","W","Amst","Bridge")
	for( i in seq_along(tmp2))
		tmp[which(df[,RegionConcentrated==tmp2[i]])]	<- i
	df[,sort1:=factor(tmp, levels=1:7, labels=tmp2) ]
	df[,sort2:=cluster.PosT]
	
	tmp					<- numeric(nrow(df))		
	tmp2				<- c("Curu","N","E","S","W","Amst")
	for( i in seq_along(tmp2))
		tmp[which(df[,RegionHospital==tmp2[i]])]	<- i
	df[,covariate:=factor(tmp, levels=1:6, labels=tmp2) ]
	xlab				<- "first pos seq"
	xlim				<- range( df[,fpos1], na.rm=1 )
	ylab				<- ""#paste("clusters BS=",thresh.bs*100," BRL.med=",thresh.brl*100," version=",gsub('/',':',signat),sep="")	
	
	#extract clusters one by one in desired order, with all the information for plotting
	tmp					<- df[,list(sort1=sort1[1], sort2=sort2[1]),by=cluster]
	clusters.sortby1	<- as.numeric( tmp[,sort1] )
	clusters.sortby2	<- as.numeric( tmp[,sort2] )
	clusters.sortby		<- order(clusters.sortby1,clusters.sortby2, decreasing=F)
	clusters			<- tmp[clusters.sortby,cluster]														
	clusters			<- lapply(clusters,function(x)	subset(df,cluster==x,c(time,covariate,rout,RegionOrigin, cluster.cex,sort1))		)
	ylim				<- c(-2,length(clusters))
	xlim[1]				<- xlim[1]-700
	
	if(0)	#plot hospital region in colors to justify ordering of plot
	{
		ncols				<- length(levels(df[,covariate]))
		cols				<- brewer.pal(ncols,"Set1")[c(1,6,3,4,5,2)]		
		cols				<- sapply(cols, function(x) my.fade.col(x,0.7))	
		
		
		file				<- paste(dir.name,paste(infile,"RegionHospital_",gsub('/',':',signat),".pdf",sep=''),sep='/')
		cat(paste("\nwrite plot to",file))
		pdf(file,width=7,height=14)
		plot(1,1,type='n',bty='n',xlim=xlim,ylim=ylim,ylab=ylab, xaxt='n',yaxt='n',xlab=xlab)
		axis.Date(1, seq.Date(xlim[1],xlim[2],by="year"), labels = TRUE )
		dummy	<- sapply(seq_along(clusters),function(i)
				{							
					cluster.cex	<- cex.points  * clusters[[i]][1,cluster.cex]
					cluster.ix	<- order(as.numeric(clusters[[i]][,time]))
					cluster.x	<- as.numeric(clusters[[i]][,time])[cluster.ix]										
					cluster.y	<- rep(i,nrow(clusters[[i]]))
					cluster.z	<- as.numeric(clusters[[i]][,covariate])[cluster.ix]
					#print(clusters[[i]][,covariate])
					#if(cluster.z[1]!=4)
					#lines( c(cluster.x[1],xlim[2]), rep(i,2), col=cols[ cluster.z[1] ] )					
					points( cluster.x, cluster.y, col=cols[ cluster.z ], pch=19, cex=cluster.cex )										
				})	
		c.lines	<- sapply(clusters,function(x) x[1,sort1])
		c.lines	<- which(diff(as.numeric(c.lines))!=0)
		lapply(c.lines,function(x)		abline(h=x+0.5, lty=3, lwd=0.75, col="grey50")			)
		
		mtext(text=levels(df[,covariate]), 2, at=c(0,c.lines[-length(c.lines)]) + diff( c(0,c.lines) )/2)
		legend("topleft",bty='n',pt.bg=cols,pch=21,legend=levels(df[,covariate]), col=rep("transparent",ncols))				
		dev.off()	
	}
	if(0)	#plot exposure group in colors
	{
		ncols				<- length(levels(df[,rout]))
		cols				<- c(brewer.pal(ncols-1,"Set1"),"grey50")		
		cols				<- sapply(cols, function(x) my.fade.col(x,0.7))[c(1,3,2,4)]	
				
		file				<- paste(dir.name,paste(infile,"rout_",gsub('/',':',signat),".pdf",sep=''),sep='/')
		cat(paste("\nwrite plot to",file))
		pdf(file,width=7,height=14)
		plot(1,1,type='n',bty='n',xlim=xlim,ylim=ylim,ylab=ylab, xaxt='n',yaxt='n',xlab=xlab)
		axis.Date(1, seq.Date(xlim[1],xlim[2],by="year"), labels = TRUE )
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
				})	
		c.lines	<- sapply(clusters,function(x) x[1,sort1])
		c.lines	<- which(diff(as.numeric(c.lines))!=0)
		lapply(c.lines,function(x)		abline(h=x+0.5, lty=3, lwd=0.75, col="grey50")			)
		
		mtext(text=levels(df[,covariate]), 2, at=c(0,c.lines[-length(c.lines)]) + diff( c(0,c.lines) )/2)
		legend("topleft",bty='n',pt.bg=cols,pch=21,legend=levels(df[,rout]), col=rep("transparent",ncols))				
		dev.off()	
	}
	if(1)	#plot region origin in colors
	{
		ncols				<- length(levels(df[,RegionOrigin]))
		cols				<- c("darkblue","firebrick1","darkorange","Cyan","grey50")	#c(brewer.pal(ncols-1,"Set1"),"grey50") 						
		cols				<- sapply(cols, function(x) my.fade.col(x,0.7))	
		
		file				<- paste(dir.name,paste(infile,"RegionOrigin_",gsub('/',':',signat),".pdf",sep=''),sep='/')
		cat(paste("\nwrite plot to",file))
		pdf(file,width=7,height=14)
		plot(1,1,type='n',bty='n',xlim=xlim,ylim=ylim,ylab=ylab, xaxt='n',yaxt='n',xlab=xlab)
		axis.Date(1, seq.Date(xlim[1],xlim[2],by="year"), labels = TRUE )
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
				})	
		c.lines	<- sapply(clusters,function(x) x[1,sort1])
		c.lines	<- which(diff(as.numeric(c.lines))!=0)
		lapply(c.lines,function(x)		abline(h=x+0.5, lty=3, lwd=0.75, col="grey50")			)
		
		mtext(text=levels(df[,covariate]), 2, at=c(0,c.lines[-length(c.lines)]) + diff( c(0,c.lines) )/2)
		legend("topleft",bty='n',pt.bg=cols,pch=21,legend=levels(df[,RegionOrigin]), col=rep("transparent",ncols))				
		dev.off()	
	}
	stop()
}

project.bezemer2013a.figs.v130821<- function()
{
	require(data.table)
	require(RColorBrewer)	 
	
	dir.name			<- "/Users/Oliver/duke/2013_HIV_NL/Bezemer2013_clusters"	
	signat				<- "130715"
	f.name				<- paste("Bezemer2014_clusters_",signat,".csv",sep='')
	df					<- read.csv(paste(dir.name,f.name,sep='/'), stringsAsFactors=FALSE)
	dir.name			<- "/Users/Oliver/duke/2013_HIV_NL/Bezemer2013_clusters"
	CLUSTER.N.THRESHOLD	<- 10
	ylab				<- ""
	yheight				<- 1.5
	cex.points			<- 0.75
	xlab				<- "first date of HIV-1 diagnosis"
	outfile				<- "BezemerFig1"
	outsignat			<- "130821"
	#deal with missing data and get into reasonable shape
	#@TODO Daniela, see if I interpret all the fields correctly, especially what is NA
	df[df[,"fpos1"]=="","fpos1"]						<- NA
	df[,"fpos1"]										<- as.Date(df[,"fpos1"], format="%d/%m/%Y")
	df[which(is.na(df[,"hcv"])),"hcv"]					<- 2
	df[,"hcv"]											<- factor(df[,"hcv"],levels=c(0,1,2),labels=c("No","Yes","unknown"))
	df[,"patient"]										<- factor(df[,"patient"])
	tmp													<- !is.na(df[,"RegionOrigin"]) & (df[,"RegionOrigin"]=="" | df[,"RegionOrigin"]=="U") 		
	df[tmp,"RegionOrigin"]								<- NA
	df[,"RegionOrigin"]									<- factor(df[,"RegionOrigin"])	
	tmp													<- unclass(df[,"RegionOrigin"]	)
	code												<- seq_along(attr(tmp,"levels"))
	names(code)											<- attr(tmp,"levels")
	print(code)
	tmp[ tmp==code[c("NL")] ]							<- 20
	tmp[ tmp==code[c("NLseAnt")] ]						<- 21
	tmp[ tmp==code[c("Lat")] ]							<- 22
	tmp[ tmp%in%code[c("EUC","EUO","EUW")] ]			<- 23
	tmp[ tmp%in%code[c("AUS","NAM","OAP","SSA","ZAz")] ]<- 24
	df[,"RegionOrigin"]									<- factor(tmp, levels=as.character(20:24), labels=c("NL","NLseAnt","Latin","Europe","other"))
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
	df	<- data.table(df, key="clustername")
	setnames(df, c("clustername","hregion"),c("cluster","RegionHospital"))
	str(df)
	#
	#	determine median fpos1 time for cluster and add PLUS determine cex for each cluster and add
	#
	tmp			<- df[,	list(	cluster.PosT=	mean(fpos1, na.rm=T), 
								cluster.cex=	log(sequencesperCluster[1]),
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
	clusters.small		<- df[, list(select= (sequencesperCluster>1 & sequencesperCluster<CLUSTER.N.THRESHOLD)[1]), by= cluster]	
	clusters.small		<- merge(subset(clusters.small, select), df, by="cluster")
	cat(paste("number seq in small clusters, n=", nrow(clusters.small)))
	df					<- subset(df, !is.na(cluster) & sequencesperCluster>=CLUSTER.N.THRESHOLD)
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
	tmp2				<- c("DU","HT","MSM","mixed")
	for( i in seq_along(tmp2))
		tmp[which(df[,ExpGrConcentrated==tmp2[i]])]	<- i
	df[,sort1:=factor(tmp, levels=1:4, labels=tmp2) ]
	df[,sort2:=max(duration,na.rm=1)-duration]
	#	
	xlim				<- range( df[,fpos1], na.rm=1 )		
	#extract clusters one by one in desired order, with all the information for plotting
	tmp					<- df[,list(sort1=sort1[1], sort2=sort2[1]),by=cluster]
	clusters.sortby1	<- as.numeric( tmp[,sort1] )
	clusters.sortby2	<- as.numeric( tmp[,sort2] )
	clusters.sortby		<- order(clusters.sortby1,clusters.sortby2, decreasing=F)
	clusters			<- tmp[clusters.sortby,cluster]														
	clusters			<- lapply(clusters,function(x)	subset(df,cluster==x,c(time, rout, RegionOrigin, RegionHospital, hcv, cluster.cex, sort1))		)
	ylim				<- c(-50,length(clusters))
	xlim[1]				<- xlim[1]-700	
	
	#
	#	plot by HCV
	#
	ncols				<- length(levels(df[,hcv]))
	cols				<- c(brewer.pal(9,"Set1"),"grey50")		
	cols				<- sapply(cols, function(x) my.fade.col(x,0.8))[c(2,8,10)]	
	file				<- paste(dir.name,paste(outfile,"_hcv_",gsub('/',':',outsignat),".pdf",sep=''),sep='/')
	cat(paste("\nwrite plot to",file))
	pdf(file,width=7,height=12)
	par(mar=c(4,4,0,0))	
	plot(1,1,type='n',bty='n',xlim=xlim,ylim=ylim,ylab=ylab, xaxt='n',yaxt='n',xlab=xlab)
	axis.Date(1, seq.Date(xlim[1],xlim[2],by="year"), labels = TRUE )
	dummy	<- sapply(seq_along(clusters),function(i)
			{							
				cluster.cex	<- cex.points  * clusters[[i]][1,cluster.cex]
				cluster.ix	<- order(as.numeric(clusters[[i]][,time]))
				cluster.x	<- as.numeric(clusters[[i]][,time])[cluster.ix]										
				cluster.y	<- rep(i,nrow(clusters[[i]]))
				cluster.z	<- as.numeric(clusters[[i]][,hcv])[cluster.ix]
				#print(clusters[[i]][,covariate])
				#if(cluster.z[1]!=4)
				#lines( c(cluster.x[1],xlim[2]), rep(i,2), col=cols[ cluster.z[1] ] )					
				points( cluster.x, cluster.y, col=cols[ cluster.z ], pch=19, cex=cluster.cex )										
			})	
	
	c.lines	<- sapply(clusters,function(x) x[1,sort1])
	c.lines	<- c( which(diff(as.numeric(c.lines))!=0), length(clusters) )
	dummy	<- lapply(c(0,c.lines),function(x)		abline(h=x+0.5, lty=3, lwd=0.75, col="grey50")			)	
	mtext(text=levels(df[,sort1]), 2, at=c(0,c.lines[-length(c.lines)]) + diff( c(0,c.lines) )/2, las=1, cex=0.7)	
	legend(xlim[1]-600,ylim[2]*0.9,bty='n',pt.bg=cols,pch=21,legend=levels(df[,hcv]), col=rep("transparent",ncols))			
	#plot small clusters	
	tmp	<- runif(nrow(clusters.small), min=-23,max=-3)
	points( as.numeric(clusters.small[,fpos1]), tmp, col=cols[ as.numeric(clusters.small[,hcv]) ], pch=19, cex= cex.points*0.4 )
	mtext(text="patients\nclusters<10", 2, at=mean(c(-23,-3)), las=1, cex=0.7)
	#plot singletons
	tmp	<- runif(nrow(singletons), min=-47,max=-27)
	points( as.numeric(singletons[,fpos1]), tmp, col=cols[ as.numeric(singletons[,hcv]) ], pch=19, cex= cex.points*0.4 )
	mtext(text="patients\nnot\nclustering", 2, at=mean(c(-47,-27)), las=1, cex=0.7)	
	dev.off()
	#
	#	plot by region origin
	#
	ncols				<- length(levels(df[,RegionOrigin]))
	cols				<- c(brewer.pal(9,"Set1"),"grey50")		
	cols				<- sapply(cols, function(x) my.fade.col(x,0.8))[c(2,5,3,7,10)]	
	file				<- paste(dir.name,paste(outfile,"_RegionOrigin_",gsub('/',':',outsignat),".pdf",sep=''),sep='/')
	cat(paste("\nwrite plot to",file))
	pdf(file,width=7,height=12)
	par(mar=c(4,4,0,0))	
	plot(1,1,type='n',bty='n',xlim=xlim,ylim=ylim,ylab=ylab, xaxt='n',yaxt='n',xlab=xlab)
	axis.Date(1, seq.Date(xlim[1],xlim[2],by="year"), labels = TRUE )
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
			})	
	
	c.lines	<- sapply(clusters,function(x) x[1,sort1])
	c.lines	<- c( which(diff(as.numeric(c.lines))!=0), length(clusters) )
	dummy	<- lapply(c(0,c.lines),function(x)		abline(h=x+0.5, lty=3, lwd=0.75, col="grey50")			)	
	mtext(text=levels(df[,sort1]), 2, at=c(0,c.lines[-length(c.lines)]) + diff( c(0,c.lines) )/2, las=1, cex=0.7)	
	legend(xlim[1]-600,ylim[2]*0.9,bty='n',pt.bg=cols,pch=21,legend=levels(df[,RegionOrigin]), col=rep("transparent",ncols))			
	#plot small clusters	
	tmp	<- runif(nrow(clusters.small), min=-23,max=-3)
	points( as.numeric(clusters.small[,fpos1]), tmp, col=cols[ as.numeric(clusters.small[,RegionOrigin]) ], pch=19, cex= cex.points*0.4 )
	mtext(text="patients\nclusters<10", 2, at=mean(c(-23,-3)), las=1, cex=0.7)
	#plot singletons
	tmp	<- runif(nrow(singletons), min=-47,max=-27)
	points( as.numeric(singletons[,fpos1]), tmp, col=cols[ as.numeric(singletons[,RegionOrigin]) ], pch=19, cex= cex.points*0.4 )
	mtext(text="patients\nnot\nclustering", 2, at=mean(c(-47,-27)), las=1, cex=0.7)	
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
	plot(1,1,type='n',bty='n',xlim=xlim,ylim=ylim,ylab=ylab, xaxt='n',yaxt='n',xlab=xlab)
	axis.Date(1, seq.Date(xlim[1],xlim[2],by="year"), labels = TRUE )
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
			})	
	
	c.lines	<- sapply(clusters,function(x) x[1,sort1])
	c.lines	<- c( which(diff(as.numeric(c.lines))!=0), length(clusters) )
	dummy	<- lapply(c(0,c.lines),function(x)		abline(h=x+0.5, lty=3, lwd=0.75, col="grey50")			)	
	mtext(text=levels(df[,sort1]), 2, at=c(0,c.lines[-length(c.lines)]) + diff( c(0,c.lines) )/2, las=1, cex=0.7)	
	legend(xlim[1]-600,ylim[2]*0.9,bty='n',pt.bg=cols,pch=21,legend=levels(df[,rout]), col=rep("transparent",ncols))			
	#plot small clusters	
	tmp	<- runif(nrow(clusters.small), min=-23,max=-3)
	points( as.numeric(clusters.small[,fpos1]), tmp, col=cols[ as.numeric(clusters.small[,rout]) ], pch=19, cex= cex.points*0.4 )
	mtext(text="patients\nclusters<10", 2, at=mean(c(-23,-3)), las=1, cex=0.7)
	#plot singletons
	tmp	<- runif(nrow(singletons), min=-47,max=-27)
	points( as.numeric(singletons[,fpos1]), tmp, col=cols[ as.numeric(singletons[,rout]) ], pch=19, cex= cex.points*0.4 )
	mtext(text="patients\nnot\nclustering", 2, at=mean(c(-47,-27)), las=1, cex=0.7)	
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
	clusters			<- lapply(clusters,function(x)	subset(df,cluster==x,c(time, rout, RegionOrigin, RegionHospital, hcv, cluster.cex, sort1))		)
	ylim				<- c(-50,length(clusters))
	xlim[1]				<- xlim[1]-700	

	ncols				<- length(levels(df[,RegionHospital]))
	cols				<- c(brewer.pal(9,"Set1"),"grey50")		
	cols				<- sapply(cols, function(x) my.fade.col(x,0.8))[c(1,3,2,4,5,6,7,10)]	
	file				<- paste(dir.name,paste(outfile,"_RegionHospital_",gsub('/',':',outsignat),".pdf",sep=''),sep='/')
	cat(paste("\nwrite plot to",file))
	pdf(file,width=7,height=12)
	par(mar=c(4,4,0,0))	
	plot(1,1,type='n',bty='n',xlim=xlim,ylim=ylim,ylab=ylab, xaxt='n',yaxt='n',xlab=xlab)
	axis.Date(1, seq.Date(xlim[1],xlim[2],by="year"), labels = TRUE )
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
			})	
	
	c.lines	<- sapply(clusters,function(x) x[1,sort1])
	c.lines	<- c( which(diff(as.numeric(c.lines))!=0), length(clusters) )
	dummy	<- lapply(c(0,c.lines),function(x)		abline(h=x+0.5, lty=3, lwd=0.75, col="grey50")			)	
	mtext(text=levels(df[,sort1]), 2, at=c(0,c.lines[-length(c.lines)]) + diff( c(0,c.lines) )/2, las=1, cex=0.5)	
	legend(xlim[1]-600,ylim[2]*1,bty='n',pt.bg=cols,pch=21,legend=levels(df[,RegionHospital]), col=rep("transparent",ncols))			
	#plot small clusters	
	tmp	<- runif(nrow(clusters.small), min=-23,max=-3)
	points( as.numeric(clusters.small[,fpos1]), tmp, col=cols[ as.numeric(clusters.small[,RegionHospital]) ], pch=19, cex= cex.points*0.4 )
	mtext(text="patients\nclusters<10", 2, at=mean(c(-23,-3)), las=1, cex=0.7)
	#plot singletons
	tmp	<- runif(nrow(singletons), min=-47,max=-27)
	points( as.numeric(singletons[,fpos1]), tmp, col=cols[ as.numeric(singletons[,RegionHospital]) ], pch=19, cex= cex.points*0.4 )
	mtext(text="patients\nnot\nclustering", 2, at=mean(c(-47,-27)), las=1, cex=0.7)	
	dev.off()
	
}

project.bezemer2013a.figs.v130528<- function()
{
	require(data.table)
	require(RColorBrewer)	 
	
	dir.name		<- "/Users/Oliver/duke/2013_HIV_NL/Bezemer2013_clusters"	
	signat			<- "130528"
	f.name			<- paste("Bezemer2014_clusters_",signat,".csv",sep='')
	df				<- read.csv(paste(dir.name,f.name,sep='/'), stringsAsFactors=FALSE)
	
	#deal with missing data and get into reasonable shape
	#@TODO Daniela, see if I interpret all the fields correctly, especially what is NA
	df[df[,"fpos1"]=="","fpos1"]				<- NA
	df[,"fpos1"]								<- as.Date(df[,"fpos1"], format="%d/%m/%Y")	
	df[,"patient"]								<- factor(df[,"patient"])
	tmp											<- !is.na(df[,"RegionOrigin"]) & (df[,"RegionOrigin"]=="" | df[,"RegionOrigin"]=="U") 		
	df[tmp,"RegionOrigin"]						<- NA
	df[,"RegionOrigin"]							<- factor(df[,"RegionOrigin"])
	df[which(df[,"Country"]==""),"Country"]		<- NA
	df[,"Country"]								<- factor(df[,"Country"])
	df[which(df[,"hregion"]==""),"hregion"]		<- "unknown"
	df[,"hregion"]								<- factor(df[,"hregion"], levels=c("CU","Ea","NH","No","So","UT","ZH","unknown"))
	df[which(df[,"sexy"]==""),"sexy"]			<- NA
	df[,"sexy"]									<- factor(df[,"sexy"])
	df[which(df[,"rout"]==""),"rout"]			<- "unknown"
	df[which(df[,"rout"]=="U"),"rout"]			<- "unknown"
	df[,"rout"]									<- factor(df[,"rout"])
	#print(levels(df[,"Country"]))
	#print(levels(df[,"hregion"]))
	#print(levels(df[,"sexy"]))
	#print(levels(df[,"rout"]))
	#print(summary(df[,"fpos1"]))
	df	<- data.table(df, key="clustername")
	str(df)
	if(0)	#evaluate distribution by "Country" and collect countries from which there is above level migration to NL 
	{
		tmp<- sort(table(df[,Country]),decreasing=TRUE)
		barplot(tmp, cex.names=.75)
		tmp<- rev(tmp)/sum(tmp)
		levels<- c(0.01,0.05)
		countries.migration.abovelevel	<- lapply(levels,function(level)
											{
												names(tmp)[cumsum(tmp)>level]
											})
		names(countries.migration.abovelevel)<- levels							
		dput(countries.migration.abovelevel)				
	}
	dir.name			<- "/Users/Oliver/duke/2013_HIV_NL/Bezemer2013_clusters"
	CLUSTER.N.THRESHOLD	<- 10
	ylab				<- "transmission network ID"
	yheight				<- 1.5
	cex.points			<- 0.75
	#networks with >= 10 sequences by date of diagnosis, ordered by duration of the network (added)
	#bottom line: singletons and smaller clusters as two summary networks on the bottom line.
	 
	#C, by country origin, blue NL, orange NLseAnt, green LAtin, black Europe. 
	#
	if(0)
	{
		#input
		f.name					<- paste("BezemerFig1_RegionOrigin_wSingletons_",signat,".pdf",sep='')
		#reduce risk group as desired and set as covariate
		tmp													<- unclass(df[,RegionOrigin]	)
		code												<- seq_along(attr(tmp,"levels"))
		names(code)											<- attr(tmp,"levels")
		print(code)
		tmp[ tmp==code[c("NL")] ]							<- 20
		tmp[ tmp==code[c("NLseAnt")] ]						<- 21
		tmp[ tmp==code[c("Lat")] ]							<- 22
		tmp[ tmp%in%code[c("EUC","EUO","EUW")] ]			<- 23
		tmp[ tmp%in%code[c("AUS","NAM","OAP","SSA","ZAz")] ]<- 24
		tmp													<- factor(tmp, levels=as.character(20:24), labels=c("NL","NLseAnt","Latin","Europe","other"))
		df[,covariate:=tmp ]
		df[,covariateNA:=tmp]
		
		
		print( attr(unclass(df[,RegionOrigin]),"levels") )
		print( attr(unclass(df[,covariate]),"levels") )
		
		xlim			<- range( df[,fpos1], na.rm=1 )
		zlim			<- range( df[,duration], na.rm=1 )
		#select singletons
		singletons		<- subset(df,clustername==0,c(fpos1,covariate,duration))
		#rm singletons from data table
		df				<- subset(df,clustername!=0,)
		#select clusters with less than CLUSTER.N.THRESHOLD non-NA fpos1 members
		cluster.n		<- df[,length(which(!is.na(fpos1))),by=clustername]		
		clusters.small	<- df[ subset(cluster.n,V1<CLUSTER.N.THRESHOLD,clustername) ]		
		#select clusters with more than CLUSTER.N.THRESHOLD non-NA fpos1 members
		df				<- df[ subset(cluster.n,V1>=CLUSTER.N.THRESHOLD,clustername) ]
		#sort clusters by duration -- this is the same as median duration
		
		clusters.sortby1<- as.vector( df[,as.numeric(any(hregion=="CU", na.rm=1)), by=clustername]$V1 )
		clusters.sortby2<- as.vector( df[,median(duration, na.rm=1), by=clustername]$V1 )
		clusters.sortby	<- order(clusters.sortby1,clusters.sortby2, decreasing=T)		
		#extract clusters one by one
		clusters		<- unique( df[,clustername] )[clusters.sortby]		
		clusters.levels	<- levels(df[,covariate])			
		clusters		<- lapply(clusters,function(x)
				{
					subset(df,clustername==x,c(fpos1,covariate,duration,hregion))
				})
		ylim			<- c(-15,length(clusters))
		cols			<- c(brewer.pal(length(clusters.levels)-1, "Set1"), "grey30")[c(2,1,3,4,5)]		#add grey for "unknown"
		#cols			<- c("green","red","blue","grey30")
		cols			<- sapply(cols, function(x) my.fade.col(x,0.7))
		pch				<- c(rep(19,length(clusters.levels)))
		xlim[1]			<- xlim[1]-700
		file			<- paste(dir.name,f.name,sep='/')
		cat(paste("\nwrite plot to",file))
		pdf(file,width=7,height=yheight*7)
		plot(1,1,type='n',bty='n',xlim=xlim,ylim=ylim,ylab=ylab, xaxt='n',xlab="")
		axis.Date(1, seq.Date(xlim[1],xlim[2],by="year"), labels = TRUE )
		#axis(3, )
		sapply(seq_along(clusters),function(i)
				{							
					cluster.cex	<- (1 - nrow(subset(clusters[[i]],hregion=="unknown")) / nrow(clusters[[i]])) * cex.points
					cluster.ix	<- order(as.numeric(clusters[[i]][,fpos1]))
					cluster.x	<- as.numeric(clusters[[i]][,fpos1])[cluster.ix]										
					cluster.y	<- rep(i,nrow(clusters[[i]]))
					cluster.z	<- as.numeric(clusters[[i]][,covariate])[cluster.ix]
					#print(clusters[[i]][,covariate])
					#if(cluster.z[1]!=4)
					lines( c(cluster.x[1],xlim[2]), rep(i,2), col=cols[ cluster.z[1] ] )					
					points( cluster.x, cluster.y, col=cols[ cluster.z ], pch=pch[1], cex=cluster.cex )										
				})		
		#plot small clusters
		abline(h=-5, lty=3,col="grey50")
		points( as.numeric(clusters.small[,fpos1]), rep(-5,nrow(clusters.small)), col=cols[ as.numeric(clusters.small[,covariate]) ], pch=pch[1], cex=0.75 )
		#plot singletons
		abline(h=-10, lty=4,col="grey50")
		points( as.numeric(singletons[,fpos1]), rep(-10,nrow(singletons)), col=cols[ as.numeric(singletons[,covariate]) ], pch=pch[1], cex=0.75 )
		
		pch[pch==19]	<- 21
		col				<- c(rep("transparent",length(clusters.levels)),"black")
		legend("topleft",bty='n',pt.bg=cols,pch=pch,legend=c(levels(df[,covariate])), col=col)
		legend("bottomleft", bty='n', lty=c(4), col="grey50", legend=c("singletons"))
		legend("bottomright", bty='n', lty=c(3), col="grey50", legend=c("networks<10"))
		dev.off()			
	}
	#B, in blue hregion 'NH', in orange in 'CU', in green in other provinces.
	if(0)	#los alamos left out, figure x: fpos1 y: clusters, coloured by hregion
	{
		#input
		f.name						<- paste("BezemerFig1_hregion_wSingletons_",signat,".pdf",sep='')
		#reduce risk group as desired and set as covariate
		tmp												<- unclass(df[,hregion]	)
		code											<- seq_along(attr(tmp,"levels"))
		names(code)										<- attr(tmp,"levels")						
		tmp[ tmp%in%code[c("Ea","No","So","UT","ZH")] ]	<- 12
		tmp[ tmp==code[c("NH")] ]						<- 10
		tmp[ tmp==code[c("CU")] ]						<- 11
		tmp[ tmp==code[c("unknown")] ]					<- 13	
		tmp												<- factor(tmp, levels=as.character(10:13), labels=c("NH","CU","other","unknown"))
		df[,covariate:=tmp ]
		tmp[ which(tmp=="unknown") ]					<- NA
		df[,covariateNA:=tmp]
		#print(df)
		
		xlim			<- range( df[,fpos1], na.rm=1 )
		zlim			<- range( df[,duration], na.rm=1 )
		#select singletons
		singletons		<- subset(df,clustername==0,c(fpos1,covariate,duration))
		#rm singletons from data table
		df				<- subset(df,clustername!=0,)
		#select clusters with less than CLUSTER.N.THRESHOLD non-NA fpos1 members
		cluster.n		<- df[,length(which(!is.na(fpos1))),by=clustername]		
		clusters.small	<- df[ subset(cluster.n,V1<CLUSTER.N.THRESHOLD,clustername) ]		
		#select clusters with more than CLUSTER.N.THRESHOLD non-NA fpos1 members
		df				<- df[ subset(cluster.n,V1>=CLUSTER.N.THRESHOLD,clustername) ]
		#sort clusters by duration -- this is the same as median duration
		
		clusters.sortby1<- as.vector( df[,as.numeric(any(covariateNA=="CU", na.rm=1)), by=clustername]$V1 )
		clusters.sortby2<- as.vector( df[,median(duration, na.rm=1), by=clustername]$V1 )
		clusters.sortby	<- order(clusters.sortby1,clusters.sortby2, decreasing=T)		
		#extract clusters one by one
		clusters		<- unique( df[,clustername] )[clusters.sortby]		
		clusters.levels	<- levels(df[,covariate])			
		clusters		<- lapply(clusters,function(x)
				{
					subset(df,clustername==x,c(fpos1,covariate,duration,hregion))
				})
		ylim			<- c(-15,length(clusters))
		cols			<- c(brewer.pal(length(clusters.levels)-1, "Set1"), "grey30")[c(2,1,3,4)]		#add grey for "unknown"
		#cols			<- c("green","red","blue","grey30")
		cols			<- sapply(cols, function(x) my.fade.col(x,0.7))
		pch				<- c(rep(19,length(clusters.levels)))
		xlim[1]			<- xlim[1]-700
		file			<- paste(dir.name,f.name,sep='/')
		cat(paste("\nwrite plot to",file))
		pdf(file,width=7,height=yheight*7)
		plot(1,1,type='n',bty='n',xlim=xlim,ylim=ylim,ylab=ylab, xaxt='n',xlab="")
		axis.Date(1, seq.Date(xlim[1],xlim[2],by="year"), labels = TRUE )
		#axis(3, )
		sapply(seq_along(clusters),function(i)
				{							
					cluster.cex	<- (1 - nrow(subset(clusters[[i]],hregion=="unknown")) / nrow(clusters[[i]])) * cex.points
					cluster.ix	<- order(as.numeric(clusters[[i]][,fpos1]))
					cluster.x	<- as.numeric(clusters[[i]][,fpos1])[cluster.ix]										
					cluster.y	<- rep(i,nrow(clusters[[i]]))
					cluster.z	<- as.numeric(clusters[[i]][,covariate])[cluster.ix]
					#if(cluster.z[1]!=4)
					lines( c(cluster.x[1],xlim[2]), rep(i,2), col=cols[ cluster.z[1] ] )					
					points( cluster.x, cluster.y, col=cols[ cluster.z ], pch=pch[1], cex=cluster.cex )										
				})		
		#plot small clusters
		abline(h=-5, lty=3,col="grey50")
		points( as.numeric(clusters.small[,fpos1]), rep(-5,nrow(clusters.small)), col=cols[ as.numeric(clusters.small[,covariate]) ], pch=pch[1], cex=0.75 )
		#plot singletons
		abline(h=-10, lty=4,col="grey50")
		points( as.numeric(singletons[,fpos1]), rep(-10,nrow(singletons)), col=cols[ as.numeric(singletons[,covariate]) ], pch=pch[1], cex=0.75 )
		
		pch[pch==19]	<- 21
		col				<- c(rep("transparent",length(clusters.levels)),"black")
		legend("topleft",bty='n',pt.bg=cols,pch=pch,legend=c(levels(df[,covariate])), col=col)
		legend("bottomleft", bty='n', lty=c(4), col="grey50", legend=c("singletons"))
		legend("bottomright", bty='n', lty=c(3), col="grey50", legend=c("networks<10"))
		dev.off()					
	}
	if(1)	#A, By riskgroup, blue MSM, red HT, green DU, black other. Annotate somehow those with high percentage / number of online retrieved seqs.
	{
		#output file
		f.name						<- paste("BezemerFig1_rout_wSingletons_",signat,".pdf",sep='')				
		#reduce risk group as desired and set as covariate
		tmp											<- unclass(df[,rout]	)
		code										<- seq_along(attr(tmp,"levels"))
		names(code)									<- attr(tmp,"levels")
		tmp[ tmp%in%code[c("BL","CH","unknown")] ]	<- code["unknown"]
		tmp											<- factor(tmp, levels=as.character(code[c("MSM","HT","DU","unknown")]), labels=c("MSM","HT","DU","other/unknown"))
		df[,covariate:=tmp ]
		tmp											<- df[,covariate]
		tmp[ which(tmp=="other/unknown") ]			<- NA
		df[,covariateNA:=tmp]
		
		xlim			<- range( df[,fpos1], na.rm=1 )
		zlim			<- range( df[,duration], na.rm=1 )
		#select singletons
		singletons		<- subset(df,clustername==0,c(fpos1,covariate,duration))
		#rm singletons from data table
		df				<- subset(df,clustername!=0,)
		#select clusters with less than CLUSTER.N.THRESHOLD non-NA fpos1 members
		cluster.n		<- df[,length(which(!is.na(fpos1))),by=clustername]		
		clusters.small	<- df[ subset(cluster.n,V1<CLUSTER.N.THRESHOLD,clustername) ]		
		#select clusters with more than CLUSTER.N.THRESHOLD non-NA fpos1 members
		df				<- df[ subset(cluster.n,V1>=CLUSTER.N.THRESHOLD,clustername) ]
		#sort clusters by duration -- this is the same as median duration
		clusters.sortby1<- as.vector( df[,as.numeric(any(hregion=="CU", na.rm=1)), by=clustername]$V1 )
		#clusters.sortby1<- as.vector( df[,median(as.numeric(covariateNA), na.rm=1), by=clustername]$V1 )
		clusters.sortby2<- as.vector( df[,median(duration, na.rm=1), by=clustername]$V1 )
		clusters.sortby	<- order(clusters.sortby1,clusters.sortby2, decreasing=T)
		
		#clusters.sortby	<- sort( as.vector( df[,median(duration, na.rm=1), by=clustername]$V1 ), index.return=T, decreasing=T)$ix
		clusters		<- unique( df[,clustername] )[clusters.sortby]
		#extract clusters one by one
		clusters.levels	<- levels(df[,covariate])			
		clusters		<- lapply(clusters,function(x)
				{
					subset(df,clustername==x,c(fpos1,covariate,duration,hregion))
				})
		ylim			<- c(-15,length(clusters))
		cols			<- c(brewer.pal(length(clusters.levels)-1, "Set1"), "grey30")[c(2,1,3,4)]		#add grey for "unknown"
		#cols			<- c("green","red","blue","grey30")
		cols			<- sapply(cols, function(x) my.fade.col(x,0.7))
		pch				<- c(rep(19,length(clusters.levels)))
		xlim[1]			<- xlim[1]-700
		file			<- paste(dir.name,f.name,sep='/')
		cat(paste("\nwrite plot to",file))
		pdf(file,width=7,height=yheight*7)
		plot(1,1,type='n',bty='n',xlim=xlim,ylim=ylim,ylab=ylab, xaxt='n',xlab="")
		axis.Date(1, seq.Date(xlim[1],xlim[2],by="year"), labels = TRUE )
		#axis(3, )
		sapply(seq_along(clusters),function(i)
				{							
					cluster.cex<- (1 - nrow(subset(clusters[[i]],hregion=="unknown")) / nrow(clusters[[i]])) * cex.points
					points( as.numeric(clusters[[i]][,fpos1]), rep(i,nrow(clusters[[i]])), col=cols[ as.numeric(clusters[[i]][,covariate]) ], pch=pch[1], cex=cluster.cex )										
				})		
		
		#plot small clusters
		abline(h=-5, lty=3,col="grey50")
		points( as.numeric(clusters.small[,fpos1]), rep(-5,nrow(clusters.small)), col=cols[ as.numeric(clusters.small[,covariate]) ], pch=pch[1], cex=0.75 )
		#plot singletons
		abline(h=-10, lty=4,col="grey50")
		points( as.numeric(singletons[,fpos1]), rep(-10,nrow(singletons)), col=cols[ as.numeric(singletons[,covariate]) ], pch=pch[1], cex=0.75 )
		
		pch[pch==19]	<- 21
		col				<- c(rep("transparent",length(clusters.levels)),"black")
		legend("topleft",bty='n',pt.bg=cols,pch=pch,legend=c(levels(df[,covariate])), col=col)
		legend("bottomleft", bty='n', lty=c(4), col="grey50", legend=c("singletons"))
		legend("bottomright", bty='n', lty=c(3), col="grey50", legend=c("networks<10"))
		dev.off()
		#save		
	}			
}

project.bezemer2013a.figs.v130517<- function()
{
	require(data.table)
	require(RColorBrewer)
	
	dir.name		<- "/Users/Oliver/duke/2013_HIV_NL/Bezemer2013_clusters"	
	signat			<- "130517"
	f.name			<- paste("Bezemer2014_clusters_",signat,".csv",sep='')
	df				<- read.csv(paste(dir.name,f.name,sep='/'), stringsAsFactors=FALSE)
	
	#deal with missing data and get into reasonable shape
	df[df[,"fpos1"]=="","fpos1"]			<- NA
	df[,"fpos1"]							<- as.Date(df[,"fpos1"], format="%d/%m/%Y")	
	df[,"patient"]							<- factor(df[,"patient"])	
	df[which(df[,"Country"]==""),"Country"]	<- NA
	df[,"Country"]							<- factor(df[,"Country"])
	df[which(df[,"hregion"]==""),"hregion"]	<- "unknown"
	df[,"hregion"]							<- factor(df[,"hregion"], levels=c("CU","Ea","NH","No","So","UT","ZH","unknown"))
	df[which(df[,"sexy"]==""),"sexy"]		<- NA
	df[,"sexy"]								<- factor(df[,"sexy"])
	df[which(df[,"rout"]==""),"rout"]		<- "unknown"
	df[which(df[,"rout"]=="U"),"rout"]		<- "unknown"
	df[,"rout"]								<- factor(df[,"rout"])
	#print(levels(df[,"Country"]))
	#print(levels(df[,"hregion"]))
	#print(levels(df[,"sexy"]))
	#print(levels(df[,"rout"]))
	#print(summary(df[,"fpos1"]))
	df	<- data.table(df, key="clustername")
	str(df)
	#I wondered if you'd like to make a figure by date of fpos1 and then A. coloured by rout, and B. coloured by hregion.				
	#The los alamos sequences do not have this information so you could colour by country and date of sequence (yearLA) instead, and in a second option leave them out.	
	
	dir.name			<- "/Users/Oliver/duke/2013_HIV_NL/Bezemer2013_clusters"
	CLUSTER.N.THRESHOLD	<- 10
	ylab				<- "transmission network ID"
	yheight				<- 1.5
	
	if(1)	#los alamos left out, figure x: fpos1 y: clusters, coloured by hregion
	{
		#input
		f.name						<- paste("BezemerFig1_hregion_noSingletons_",signat,".pdf",sep='')
		df[,covariate:=hregion]		#set which covariate to use
		tmp							<- df[,covariate]
		tmp[ which(tmp=="unknown") ]<- NA
		df[,covariateNA:=tmp]
		
		#create plot
		xlim			<- range( df[,fpos1], na.rm=1 )
		cluster.n		<- df[,length(fpos1),by=clustername]		
		#select clusters with more than CLUSTER.N.THRESHOLD members
		df				<- df[ subset(cluster.n,V1>=CLUSTER.N.THRESHOLD,clustername) ]
		#sort clusters by median covariate and then by median fpos	
		clusters.sortby1<- as.vector( df[,median(as.numeric(covariateNA), na.rm=1), by=clustername]$V1 )
		clusters.sortby2<- as.vector( df[,median(fpos1, na.rm=1), by=clustername]$V1 )
		clusters.sortby	<- order(clusters.sortby1,clusters.sortby2)
		clusters		<- unique( df[,clustername] )[clusters.sortby]				
		clusters.levels	<- levels(df[,covariate])			
		clusters		<- lapply(clusters,function(x)
				{
					subset(df,clustername==x,c(fpos1,covariate))
				})
		ylim			<- c(0,length(clusters))
		cols			<- c(brewer.pal(length(clusters.levels)-1, "Set3"), "grey30")		#add grey for "unknown"
		cols			<- sapply(cols, function(x) my.fade.col(x,0.7))
		pch				<- c(rep(19,length(clusters.levels)),3)		
		xlim[1]			<- xlim[1]-700
		file			<- paste(dir.name,f.name,sep='/')
		cat(paste("\nwrite plot to",file))
		pdf(file,width=7,height=yheight*7)
		plot(1,1,type='n',bty='n',xlim=xlim,ylim=ylim,xaxt='n',ylab=ylab, xlab="")
		axis.Date(1, seq.Date(xlim[1],xlim[2],by="year"), labels = TRUE )
		sapply(seq_along(clusters),function(i)
				{							
					points( as.numeric(clusters[[i]][,fpos1]), rep(i,nrow(clusters[[i]])), col=cols[ as.numeric(clusters[[i]][,covariate]) ], pch=pch[1], cex=0.75 )					
					points( median(as.numeric(clusters[[i]][,fpos1]), na.rm=1), i, col="black", pch=pch[length(pch)], cex=0.5 )
				})
		pch[pch==19]	<- 21
		col				<- c(rep("transparent",length(clusters.levels)),"black")
		legend("topleft",bty='n',pt.bg=c(cols,"black"),pch=pch,legend=c(levels(df[,covariate]),"median 1st HIV+"), col=col)
		dev.off()		
	}
	if(1)	#los alamos left out, figure x: fpos1 y: clusters, coloured by rout
	{
		#input
		f.name						<- paste("BezemerFig1_rout_noSingletons_",signat,".pdf",sep='')
		df[,covariate:=rout]		#set which covariate to use
		
		#create plot
		xlim			<- range( df[,fpos1], na.rm=1 )
		cluster.n		<- df[,length(fpos1),by=clustername]		
		#select clusters with more than CLUSTER.N.THRESHOLD members
		df				<- df[ subset(cluster.n,V1>=CLUSTER.N.THRESHOLD,clustername) ]
		#sort clusters by median fpos	
		clusters.sortby	<- sort( as.vector( df[,median(fpos1, na.rm=1), by=clustername]$V1 ), index.return=T)$ix		
		clusters		<- unique( df[,clustername] )[clusters.sortby]				
		clusters.levels	<- levels(df[,covariate])			
		clusters		<- lapply(clusters,function(x)
				{
					subset(df,clustername==x,c(fpos1,covariate))
				})
		ylim			<- c(0,length(clusters))
		cols			<- c(brewer.pal(length(clusters.levels)-1, "Set3"), "grey30")		#add grey for "unknown"
		cols			<- sapply(cols, function(x) my.fade.col(x,0.7))
		pch				<- c(rep(19,length(clusters.levels)),3)		
		xlim[1]			<- xlim[1]-700
		file			<- paste(dir.name,f.name,sep='/')
		cat(paste("\nwrite plot to",file))
		pdf(file,width=7,height=yheight*7)
		plot(1,1,type='n',bty='n',xlim=xlim,ylim=ylim,xaxt='n',ylab=ylab, xlab="")
		axis.Date(1, seq.Date(xlim[1],xlim[2],by="year"), labels = TRUE )
		sapply(seq_along(clusters)[-1],function(i)
				{							
					points( as.numeric(clusters[[i]][,fpos1]), rep(i,nrow(clusters[[i]])), col=cols[ as.numeric(clusters[[i]][,covariate]) ], pch=pch[1], cex=0.75 )					
					points( median(as.numeric(clusters[[i]][,fpos1]), na.rm=1), i, col="black", pch=pch[length(pch)], cex=0.5 )
				})
		pch[pch==19]	<- 21
		col				<- c(rep("transparent",length(clusters.levels)),"black")
		legend("topleft",bty='n',pt.bg=c(cols,"black"),pch=pch,legend=c(levels(df[,covariate]),"median 1st HIV+"), col=col)
		dev.off()
		#save		
	}	
}

project.bezemer2013a.figs<- function()
{
	#project.bezemer2013a.figs.v130517()
	#project.bezemer2013a.figs.v130528()
	#project.bezemer2013a.figs.v130714()
	#project.bezemer2013a.figs.v130715()
	project.bezemer2013a.figs.v130821()
}

project.bezemer2013b.rates.v130813<- function()
{
	require(data.table)
	require(ape)
	require(RColorBrewer)	 
	verbose			<- 1
	indir			<- "/Users/Oliver/duke/2013_HIV_NL/Bezemer2013_clusters"	
	insignat		<- "130813"
	infile			<- "BezemerRates_Timedtree6CUnetworks"
	f.name			<- paste(indir,'/',infile,'_',insignat,".nex",sep='')
	
	if(verbose) cat(paste("\nload file",f.name))
	ph	<- read.nexus(f.name)
	#
	# generate data.table from tip.labels
	#
	ph.tl			<- t(sapply(strsplit(ph$tip.label,'/'),function(x)
						{
							if(length(x)==16)	return(x)
							else if(length(x)==4)	return(  c(x[1],x[2],rep(NA,12),x[3],x[4]))
							else stop("unexpected line in file")
						}))
	colnames(ph.tl)	<- c("cluster","Patient","Trm","Sex","HCV","RegionDiagnosed","Infected","CountryInfected","XX","CountryBorn","Pos_T1","Age_T1","Age_now","Therapy","Country","PosSeq_T1")
	ph.tl			<- as.data.table(ph.tl)
	set(ph.tl,NULL,"cluster",as.numeric(as.character(ph.tl[,cluster])))
	set(ph.tl,NULL,"PosSeq_T1",as.numeric(as.character(ph.tl[,PosSeq_T1])))
	set(ph.tl,NULL,"Pos_T1",as.numeric(as.character(ph.tl[,Pos_T1])))
	set(ph.tl,NULL,"Age_T1",as.numeric(as.character(ph.tl[,Age_T1])))
	set(ph.tl,NULL,"Age_now",as.numeric(as.character(ph.tl[,Age_now])))
	# Trm
	if(verbose) cat(paste("\nsetting Trm %in% BL CH DU to NA"))
	tmp									<- as.character( ph.tl[,Trm] )
	tmp[tmp%in%c("BL","CH","DU","U")]	<- NA
	set(ph.tl,NULL,"Trm",factor(tmp, levels=c("MSM","HT"),labels=c("MSM","HET")))
	# CountryBorn
	tmp									<- as.character( ph.tl[,CountryBorn] )
	if(verbose) cat(paste("\nsetting CountryBorn %in% AN  DO  HT SR TT  VE to C(arribean)"))
	tmp[tmp%in%c("AN","CW","DO","HT","SR","TT","VE")]	<- "C"
	if(verbose) cat(paste("\nsetting CountryBorn %in% ER  GE  GR LK RO to OTH"))
	tmp[tmp%in%c("ER","GE","GR","LK","RO","")]				<- NA
	set(ph.tl, NULL, "CountryBorn", factor(tmp, levels=c("C","NL"),labels=c("C","NL")))
	# CountryDiagnosed
	tmp									<- as.character( ph.tl[,RegionDiagnosed] )
	if(verbose) cat(paste("\nsetting RegionDiagnosed %in% DR FL GL GR LB NB NH OV UT ZH to NL"))
	tmp[tmp%in%c("DR","FL","GL","GR","LB","NB","NH","OV","UT","ZH")]	<- "NL"
	tmp[tmp%in%c("CU")]	<- "C"
	ph.tl[,"CountryDiagnosed":=factor(tmp, levels=c("C","NL"),labels=c("C","NL"))]
	# CountryInfected
	tmp									<- as.character( ph.tl[,CountryInfected] )
	if(verbose) cat(paste("\nsetting CountryInfected %in% AN CW DO  HN to C"))
	tmp[tmp%in%c("AN","CW","DO","HN")]	<- "C"
	if(verbose) cat(paste("\nsetting CountryInfected %in% US GE LK ID to NA"))
	tmp[tmp%in%c("US","GE","LK","ID")]	<- NA
	set(ph.tl, NULL, "CountryInfected", factor(tmp, levels=c("C","NL"),labels=c("C","NL")))
	# MigState.null
	ph.tl[, "MigState.null":=as.character(ph.tl[,CountryDiagnosed])]
	# MigState.bridge
	ph.tl[, "MigState.bridge":=as.character(ph.tl[,CountryDiagnosed])]	
	tmp									<- c( 	which( ph.tl[, !is.na(CountryDiagnosed) & !is.na(CountryInfected) & CountryDiagnosed!=CountryInfected] ),
												which( ph.tl[, !is.na(CountryDiagnosed) & !is.na(CountryBorn) & CountryDiagnosed!=CountryBorn] ) )
	if(verbose) cat(paste("\nfound potentially bridging individuals, n=",length(tmp)))
	set(ph.tl, tmp, "MigState.bridge", "B")
	# TrmState
	ph.tl[,TrmState:= as.character(ph.tl[,Sex])]
	set(ph.tl, which(ph.tl[,TrmState=="F"]), "TrmState", "f")
	set(ph.tl, which(ph.tl[,TrmState=="M" & Trm=="HET"]), "TrmState", "m")
	set(ph.tl, which(ph.tl[,TrmState=="M" & Trm=="MSM"]), "TrmState", "msm")
	set(ph.tl, NULL, "TrmState", factor(ph.tl[,TrmState], levels=c("f","m","msm"), labels=c("f","m","msm")))	
	# TraitState.null
	tmp									<- apply( cbind( as.character( ph.tl[,MigState.null]),as.character( ph.tl[,TrmState]) ), 1, function(x) paste(x,collapse='',sep='') )
	ph.tl[,"TraitState.null":= tmp]
	# TraitState.bridge
	tmp									<- apply( cbind( as.character( ph.tl[,MigState.bridge]),as.character( ph.tl[,TrmState]) ), 1, function(x) paste(x,collapse='',sep='') )
	ph.tl[,"TraitState.bridge":= tmp]
	#print(subset(ph.tl,select=c(Patient, MigState.null, MigState.bridge, TrmState, TraitState.null, TraitState.bridge )), nrow=170)
	
	#
	# reset tip.labels
	#
	ph$tip.label						<- ph.tl[, paste(Patient,TraitState,PosSeq_T1,sep='_')]
	ph<- ladderize(ph)
	plot(ph)
	
	
}

project.bezemer2013b.rates<- function()
{
	project.bezemer2013b.rates.v130813()	
}