# Publication figures from analysis.py outputs.
suppressPackageStartupMessages({library(ggplot2);library(grid);library(jsonlite)})
args <- commandArgs(FALSE)
script <- normalizePath(sub('^--file=', '', args[grepl('^--file=', args)]), winslash='/')
root <- file.path(dirname(script), 'results')

draw_associations <- function(root) {
# Manuscript artwork from the frozen analysis; R-only graphics and export.
# Run with Rscript from any working directory. No database access or model selection.
suppressPackageStartupMessages({library(ggplot2);library(grid);library(jsonlite)})
input <- file.path(root,'primary_analysis'); out <- file.path(root,'figures/associations')
for (d in c('exports','source_data','qc')) dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
R <- fromJSON(file.path(input,'primary_results.json'))
A <- fromJSON(file.path(input,'cohort_audit.json'))
D <- fromJSON(file.path(input,'final_analysis_dataset_local.json'))
res <- fromJSON(file.path(input,'residuals_local.json'))
models <- R$models
getmodel <- function(k) models[match(k,models$label),,drop=FALSE]
stopifnot(nrow(D)==760,nrow(res)==760,A$final_n==760,A$same_subjects_values_times_as_frozen760)
stopifnot(abs(coef(lm(log(CysC)~log(FT4),D))[2]-getmodel('crude')$beta)<1e-8)
stopifnot(abs(coef(lm(y_logCysC_residual~x_logFT4_residual,res))[2]-getmodel('final_primary')$beta)<1e-8)
stopifnot(2662-6==2656,2656-1007==1649,1649-882==767,767-1-3-3==760)
file.copy(file.path(input,'primary_results.json'),file.path(out,'source_data','frozen_results.json'),overwrite=TRUE)
file.copy(file.path(input,'cohort_audit.json'),file.path(out,'source_data','cohort_counts.json'),overwrite=TRUE)
write.csv(models,file.path(out,'source_data','model_estimates.csv'),row.names=FALSE,na='')
# Observation-level values remain local and are excluded from the delivery archive.
write.csv(data.frame(log_FT4=log(D$FT4),log_CysC=log(D$CysC),res),file.path(out,'source_data','scatter_values_LOCAL_ONLY.csv'),row.names=FALSE)
blue <- '#245D72'; ink <- '#20272B'; gray <- '#60696F'; pale <- '#EDF3F5'
font <- 'Arial'; body <- 7
theme_set(theme_classic(base_size=body,base_family=font)+theme(
 axis.line=element_line(linewidth=.23,colour=ink),axis.ticks=element_line(linewidth=.23,colour=ink),
 axis.text=element_text(size=6.5,colour=ink),axis.title=element_text(size=7,colour=ink),
 plot.title=element_text(size=7,face='bold',hjust=0),plot.subtitle=element_text(size=6.5),
 plot.margin=margin(4,5,3,4),panel.grid=element_blank()))

canvas <- function(w,h){grid.newpage();pushViewport(viewport(xscale=c(0,w),yscale=c(h,0)));}
txt <- function(s,x,y,size=7,bold=FALSE,col=ink,just='left') grid.text(s,x,y,default.units='native',just=just,gp=gpar(fontfamily=font,fontsize=size,fontface=if(bold)'bold' else 'plain',col=col,lineheight=1.12))
line <- function(x1,y1,x2,y2,col=ink,lwd=.7,lty=1,arrow=FALSE) grid.lines(c(x1,x2),c(y1,y2),default.units='native',gp=gpar(col=col,fill=col,lwd=lwd,lty=lty),arrow=if(arrow)grid::arrow(length=unit(1.35,'mm'),type='closed')else NULL)
rect <- function(x,y,w,h,fill='white',col=ink,lwd=.7) grid.rect(x=unit(x,'native'),y=unit(y+h/2,'native'),width=unit(w,'mm'),height=unit(h,'mm'),just=c('left','centre'),gp=gpar(fill=fill,col=col,lwd=lwd))
box <- function(x,y,w,h,s,final=FALSE){
 rect(x,y,w,h,if(final)pale else 'white',if(final)blue else ink)
 # Check text dimensions on the actual export device, in final millimetres.
 tg <- textGrob(s,gp=gpar(fontfamily=font,fontsize=7,lineheight=1.12,fontface=if(final)'bold'else'plain'))
 stopifnot(convertWidth(grobWidth(tg),'mm',valueOnly=TRUE)<w-3,
           convertHeight(grobHeight(tg),'mm',valueOnly=TRUE)<h-2)
 txt(s,x+w/2,y+h/2,7,final,just='centre')
}
scatter <- function(x,y,title,xlab,ylab,key){
 m<-getmodel(key);d<-data.frame(x,y);itc<-mean(y)-m$beta*mean(x)
 r2<-cor(x,y)^2
 expected<-if(key=='crude')m$r2 else R$diagnostics$added_R2/(1-m$r2+R$diagnostics$added_R2)
 stopifnot(abs(r2-expected)<1e-10)
 stat<-sprintf('%s = %.3f\nHC3 P %s',if(key=='crude')'R²' else 'Partial R²',r2,
               if(m$p<.001)'< 0.001' else sprintf('= %.3f',m$p))
 ggplot(d,aes(x,y))+geom_point(size=.65,alpha=.30,colour=gray,stroke=0)+
 geom_abline(intercept=itc,slope=m$beta,colour=if(key=='crude')ink else blue,linewidth=.32)+
 annotate('text',x=-Inf,y=Inf,label=stat,hjust=-.08,vjust=1.15,size=2.35,
          family=font,colour=ink,lineheight=1.05)+
 labs(title=title,subtitle=sprintf('760 patients; slope = %.3f',m$beta),x=xlab,y=ylab)+
 scale_x_continuous(expand=expansion(mult=.06))+scale_y_continuous(expand=expansion(mult=c(.06,.20)))
}
fitlabel <- function(m) sprintf('%+.2f (%+.2f, %+.2f)',m$percent,m$low,m$high)
forestrow <- function(m,label,y,lo,hi,x0,x1,labelx=4,nx=NULL,resultx=179,primary=FALSE){
 fx<-function(z)x0+(z-lo)/(hi-lo)*(x1-x0)
 stopifnot(m$low>=lo,m$high<=hi)
 color<-if(primary)blue else ink
 txt(label,labelx,y,6.7,primary)
 if(!is.null(nx))txt(format(m$n,big.mark=','),nx,y,6.7,just='right')
 line(fx(m$low),y,fx(m$high),y,col=color,lwd=.85)
 line(fx(m$low),y-.65,fx(m$low),y+.65,col=color,lwd=.7)
 line(fx(m$high),y-.65,fx(m$high),y+.65,col=color,lwd=.7)
 grid.points(fx(m$percent),y,default.units='native',pch=if(primary)15 else 16,size=unit(if(primary)1.65 else 1.4,'mm'),gp=gpar(col=color))
 txt(fitlabel(m),resultx,y,6.6,primary,just='right')
}
axisforest <- function(y,lo,hi,x0,x1,ticks){
 fx<-function(z)x0+(z-lo)/(hi-lo)*(x1-x0)
 line(x0,y,x1,y,lwd=.65)
 for(t in ticks){line(fx(t),y,fx(t),y+1.1,lwd=.65);txt(as.character(t),fx(t),y+3,6.5,just='centre')}
}
direction <- function(){
 canvas(183,158)
 pa<-scatter(log(D$FT4),log(D$CysC),'Unadjusted association','ln FT4 (pmol/L)','ln cystatin C (mg/L)','crude')
 pb<-scatter(res$x_logFT4_residual,res$y_logCysC_residual,'Conditional association','Residual ln FT4','Residual ln cystatin C','final_primary')
 print(pa,vp=viewport(x=unit(3,'mm'),y=unit(155,'mm'),width=unit(87,'mm'),height=unit(70,'mm'),just=c('left','top')),newpage=FALSE)
 print(pb,vp=viewport(x=unit(95,'mm'),y=unit(155,'mm'),width=unit(86,'mm'),height=unit(70,'mm'),just=c('left','top')),newpage=FALSE)
 txt('a',2,2,8,TRUE);txt('b',94,2,8,TRUE)
 txt('c',2,79,8,TRUE);txt('Sequential adjustment in the same 760 patients',8,79,7,TRUE)
 txt('Model',4,89,6.7,TRUE);txt('Difference, % (95% CI)',179,89,6.7,TRUE,just='right')
 lo<--7;hi<-10;x0<-77;x1<-132
 line(x0+(0-lo)/(hi-lo)*(x1-x0),94,x0+(0-lo)/(hi-lo)*(x1-x0),139,col='#959CA0',lwd=.65,lty=2)
 keys<-c('crude','maturity','maturity_lag','plus_creatinine','final_primary')
 labels<-c('Unadjusted','Age, sex and recorded prematurity','Add signed pairing lag','Add creatinine','Add recent vasoactive prescription')
 for(i in seq_along(keys))forestrow(getmodel(keys[i]),labels[i],97+(i-1)*9.4,lo,hi,x0,x1,primary=i==5)
 axisforest(142,lo,hi,x0,x1,c(-5,0,5,10))
 txt('Cystatin C difference (%) for the fixed FT4 interquartile contrast',91.5,153,7,just='centre')
 popViewport()
}
exports<-data.frame(figure='Figure2_direction_reversal',width_mm=183,height_mm=158)
drawers<-list(direction)
for(i in 1:nrow(exports)){
 name<-exports$figure[i];w<-exports$width_mm[i]/25.4;h<-exports$height_mm[i]/25.4
 path<-file.path(out,'exports',name);draw<-drawers[[i]]
 svglite::svglite(paste0(path,'.svg'),width=w,height=h,bg='white');draw();dev.off()
 grDevices::cairo_pdf(paste0(path,'.pdf'),width=w,height=h,family=font,bg='white',onefile=FALSE);draw();dev.off()
 for(dpi in c(600,1200)){
  ragg::agg_tiff(paste0(path,'_',dpi,'dpi.tiff'),width=w,height=h,units='in',res=dpi,compression='lzw',background='white');draw();dev.off()
 }
 ragg::agg_png(paste0(path,'_preview.png'),width=w,height=h,units='in',res=220,background='white');draw();dev.off()
 # Render the actual exported PDF through R for visual quality assurance.
 invisible(pdftools::pdf_convert(paste0(path,'.pdf'),format='png',dpi=220,filenames=file.path(out,'qc',paste0(name,'_PDF_page%d.%s')),verbose=FALSE))
}
write.csv(exports,file.path(out,'qc','export_dimensions.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'qc','R_sessionInfo.txt'))
cat('All R figures rendered from frozen results.\n')

metrics<-list(n=760,unadjusted_R2=getmodel('crude')$r2,
 partial_R2=cor(res$x_logFT4_residual,res$y_logCysC_residual)^2,
 total_model_R2=getmodel('final_primary')$r2,incremental_R2=R$diagnostics$added_R2,
 unadjusted_P_HC3_t=getmodel('crude')$p,adjusted_P_HC3_t=getmodel('final_primary')$p,
 point_alpha=.30,display_P_threshold=.001,
 interpretation='Panel b shows partial R-squared, not the full-model or incremental R-squared; P values use the original HC3/t regression inference.')
write_json(metrics,file.path(out,'source_data','Figure2_annotation_statistics.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)

}

draw_sensitivity <- function(root) {
# R-only publication graphics. Run from the project root with authorized local results.
suppressPackageStartupMessages({library(ggplot2);library(grid);library(jsonlite)})
out<-file.path(root,'supplementary_analysis')
for(d in c('figures','qa'))dir.create(file.path(out,d),showWarnings=FALSE)
R<-fromJSON(file.path(root,'primary_analysis/primary_results.json'));models<-R$models
T<-fromJSON(file.path(out,'evidence/TSH_selection_audit.json'));L<-fromJSON(file.path(out,'evidence/longitudinal_diagnostics.json'),simplifyVector=FALSE)
getmodel<-function(k)models[match(k,models$label),,drop=FALSE]
font<-'Arial';ink<-'#20272B';blue<-'#245D72'
txt<-function(s,x,y,size=8,bold=FALSE,just='left')grid.text(s,x,y,default.units='native',just=just,gp=gpar(fontfamily=font,fontsize=size,fontface=if(bold)'bold'else'plain',col=ink,lineheight=1.1))
lin<-function(x1,y1,x2,y2,col=ink,lty=1)grid.lines(c(x1,x2),c(y1,y2),default.units='native',gp=gpar(col=col,lwd=.7,lty=lty))
mainkeys<-c('final_primary','nicu_interval','pairing_1h','pairing_24h','age7d','no_TSH','albumin_lactate_extended','recent_weight_extended','admission_weight_extended','severity_extended','within_person')
mainlabels<-c('Primary model','Both panels in recorded ICU interval','Reselected pairs within 1 h','Reselected pairs within 24 h','Both panels at age 7 d or later','TSH not required','Add albumin and lactate','Add recent weight','Add admission weight and age','Add physiology and lactate','Within-patient association')
allkeys<-c('historical_primary','broad_drug_definition','no_uncertain_drug_duration','no_broad_vaso_steroid','nonlinear_creatinine','albumin_base','albumin_extended','lactate_base','lactate_extended','albumin_lactate_base','albumin_lactate_extended','recent_weight_base','recent_weight_extended','admission_weight_base','admission_weight_extended','severity_base','severity_extended','severity_albumin_base','severity_albumin_extended')
alllabels<-c('Historical primary definition','Broader prescription definition','Exclude uncertain Rx duration','Exclude broad vasoactive / steroid Rx','Nonlinear creatinine','Albumin subset: base','Albumin subset: add albumin','Lactate subset: base','Lactate subset: add lactate','Joint subset: base','Joint subset: add albumin / lactate','Recent-weight subset: base','Recent-weight subset: add weight','Admission-weight subset: base + age','Admission-weight subset: add weight','Physiology subset: base','Physiology subset: add covariates','Physiology / albumin subset: base','Physiology / albumin: add covariates')
rowdata<-function(k){
 if(k=='no_TSH')return(as.data.frame(T$sensitivity[c('n','percent','low','high')]))
 if(k=='within_person'){v<-L$optimizers[[2]]$terms$within_FT4;return(data.frame(n=77,percent=v$percent,low=v$low,high=v$high))}
 getmodel(k)
}
forest<-function(keys,labels,height){
 grid.newpage();pushViewport(viewport(xscale=c(0,183),yscale=c(height,0)))
 txt('Analysis',3,6,8,TRUE);txt('n',80,6,8,TRUE,just='right');txt('Difference, % (95% CI)',180,6,8,TRUE,just='right')
 x0<-86;x1<-132;lo<--4;hi<-if(length(keys)>12)15 else 12;fx<-function(x)x0+(x-lo)/(hi-lo)*(x1-x0)
 y<-16+seq_along(keys)*8.4-8.4
 lin(fx(0),12,fx(0),max(y)+4,col='#8B9296',lty=2)
 for(i in seq_along(keys)){
  m<-rowdata(keys[i]);stopifnot(m$low>=lo,m$high<=hi)
  labelg<-textGrob(labels[i],gp=gpar(fontfamily=font,fontsize=8));stopifnot(convertWidth(grobWidth(labelg),'mm',valueOnly=TRUE)<73*170/183)
  txt(labels[i],3,y[i],8,keys[i]=='final_primary');txt(format(m$n,big.mark=','),80,y[i],8,just='right')
  col<-if(keys[i]=='final_primary')blue else ink
  lin(fx(m$low),y[i],fx(m$high),y[i],col);lin(fx(m$low),y[i]-.7,fx(m$low),y[i]+.7,col);lin(fx(m$high),y[i]-.7,fx(m$high),y[i]+.7,col)
  grid.points(fx(m$percent),y[i],default.units='native',pch=if(keys[i]=='final_primary')15 else 16,size=unit(1.6,'mm'),gp=gpar(col=col))
  txt(sprintf('%+.2f (%+.2f, %+.2f)',m$percent,m$low,m$high),180,y[i],7.8,just='right')
 }
 ay<-max(y)+8;lin(x0,ay,x1,ay)
 for(v in c(-3,0,3,6,9,12)){lin(fx(v),ay,fx(v),ay+1);txt(v,fx(v),ay+4,7.5,just='centre')}
 txt('Cystatin C difference for the fixed FT4 interquartile contrast (%)',91.5,height-4,8,just='centre');popViewport()
}
curve<-read.csv(file.path(out,'evidence/FT4_spline_curve.csv'));hist<-read.csv(file.path(out,'evidence/FT4_histogram.csv'))
theme_set(theme_classic(base_size=9,base_family='Arial')+theme(axis.text=element_text(colour=ink,size=8),axis.line=element_line(linewidth=.3),plot.margin=margin(6,10,4,6)))
p<-ggplot(curve,aes(FT4,percent))+geom_ribbon(aes(ymin=low,ymax=high),fill='#DCE9EE')+geom_hline(yintercept=0,linetype=3,colour='gray50',linewidth=.3)+geom_line(aes(linetype='Natural spline'),colour=blue,linewidth=.5)+geom_line(aes(y=linear_percent,linetype='Primary log-linear'),colour=ink,linewidth=.4)+scale_linetype_manual(values=c('Natural spline'=1,'Primary log-linear'=2))+labs(x='FT4 (pmol/L)',y='Cystatin C difference from median FT4 (%)',linetype=NULL)+theme(legend.position='top',legend.text=element_text(size=8))
h<-ggplot(hist,aes(xmin=lower,xmax=upper,ymin=0,ymax=count))+geom_rect(fill='#778990',colour='white',linewidth=.2)+labs(x='FT4 (pmol/L)',y='Patients (n)')+theme(axis.title=element_text(size=8))
spline<-function(){grid.newpage();print(p,vp=viewport(x=.5,y=.65,width=1,height=.7),newpage=FALSE);print(h,vp=viewport(x=.5,y=.15,width=1,height=.3),newpage=FALSE)}
spec<-list(list(name='Figure3_supporting_analyses',height=126,draw=function()forest(mainkeys,mainlabels,126)),list(name='FigureS1_FT4_functional_form',height=145,draw=spline),list(name='FigureS2_extended_models',height=194,draw=function()forest(allkeys,alllabels,194)))
for(s in spec){
 path<-file.path(out,'figures',s$name);w<-170/25.4;ht<-s$height*170/183/25.4
 svglite::svglite(paste0(path,'.svg'),width=w,height=ht);s$draw();dev.off()
 cairo_pdf(paste0(path,'.pdf'),width=w,height=ht,family=font);s$draw();dev.off()
 ragg::agg_tiff(paste0(path,'_600dpi.tiff'),width=w,height=ht,units='in',res=600,compression='lzw',background='white');s$draw();dev.off()
 ragg::agg_png(paste0(path,'_preview.png'),width=w,height=ht,units='in',res=220,background='white');s$draw();dev.off()
 invisible(pdftools::pdf_convert(paste0(path,'.pdf'),format='png',dpi=180,filenames=file.path(out,'qa',paste0(s$name,'_PDF_page%d.%s')),verbose=FALSE))
}
writeLines(capture.output(sessionInfo()),file.path(out,'evidence/R_sessionInfo.txt'))
write.csv(do.call(rbind,lapply(seq_along(mainkeys),function(i){z<-rowdata(mainkeys[i]);data.frame(label=mainlabels[i],n=z$n,percent=z$percent,low=z$low,high=z$high)})),file.path(out,'evidence/Figure3_source.csv'),row.names=FALSE)
cat('R figure export complete\n')

}

draw_dag <- function(root) {
# Hypothesis-only supplementary DAG. R is the sole graphics/export backend.
suppressPackageStartupMessages({library(grid);library(jsonlite)})
out<-file.path(root,'selection_analysis')
dir.create(file.path(out,'figures'),showWarnings=FALSE);dir.create(file.path(out,'qa'),showWarnings=FALSE)
nodes<-data.frame(id=c('M','FT4','GFR','Cr','CysC'),x=c(85,25,145,85,145),y=c(19,47,47,77,99),w=c(52,33,39,38,36),h=c(15,13,15,14,13))
edges<-data.frame(from=c('M','M','M','GFR','GFR','FT4'),to=c('FT4','GFR','Cr','Cr','CysC','CysC'),hypothesis=c(FALSE,FALSE,FALSE,FALSE,FALSE,TRUE))
# Verify directed acyclicity and that the depicted backdoor path has Cr as a collider.
walk<-function(node,trail=character()){stopifnot(!node %in% trail);for(child in edges$to[edges$from==node])walk(child,c(trail,node))}
for(n in nodes$id)walk(n)
stopifnot(all(c('M>FT4','M>Cr','GFR>Cr','GFR>CysC') %in% paste(edges$from,edges$to,sep='>')))
ink<-'#20272B';blue<-'#245D72';amber<-'#8B5D2E';font<-'Arial'
txt<-function(s,x,y,size=8,bold=FALSE,col=ink,just='centre'){grid.text(s,x,y,default.units='native',just=just,gp=gpar(fontfamily=font,fontsize=size,fontface=if(bold)'bold'else'plain',col=col,lineheight=1.12))}
arrowline<-function(x,y,col=ink,lty=1){grid.lines(x,y,default.units='native',arrow=arrow(type='closed',length=unit(1.6,'mm')),gp=gpar(col=col,fill=col,lwd=.8,lty=lty,linejoin='round'))}
box<-function(x,y,w,h,label,unmeasured=FALSE,adjusted=FALSE){
 grid.rect(x,y,w,h,default.units='native',gp=gpar(fill=if(adjusted)'#F6EFE5'else'white',col=if(adjusted)amber else ink,lwd=.8,lty=if(unmeasured)2 else 1))
 if(adjusted)grid.rect(x,y,w-1.5,h-1.5,default.units='native',gp=gpar(fill=NA,col=amber,lwd=.6))
 tg<-textGrob(label,gp=gpar(fontfamily=font,fontsize=8,lineheight=1.12));stopifnot(convertWidth(grobWidth(tg),'mm',TRUE)<w-3,convertHeight(grobHeight(tg),'mm',TRUE)<h-2)
 txt(label,x,y,8)
}
draw<-function(){
 grid.newpage();pushViewport(viewport(xscale=c(0,170),yscale=c(170,0)))
 txt('a  Hypothesized relations',4,5,9,TRUE,just='left')
 # Endpoints terminate at node boundaries, without crossing label boxes.
 arrowline(c(69,38),c(26.5,40.5));arrowline(c(101,131),c(26.5,39.5));arrowline(c(85,85),c(26.5,70))
 arrowline(c(133,104),c(54.5,70));arrowline(c(145,145),c(54.5,92.5))
 arrowline(c(25,25,127),c(53.5,99,99),col=blue,lty=2)
 box(85,19,52,15,'Gestational / developmental\nmaturation (M)',unmeasured=TRUE)
 box(25,47,33,13,'FT4')
 box(145,47,39,15,'True GFR',unmeasured=TRUE)
 box(85,77,38,14,'Creatinine\n(conditioned on)',adjusted=TRUE)
 box(145,99,36,13,'Cystatin C')
 txt('Possible non-filtration pathway',70,105,7.5,col=blue)
 txt('b  One noncausal backdoor path in panel a',4,119,9,TRUE,just='left')
 xx<-c(17,51,85,119,153);ww<-c(25,26,31,29,25);yy<-137
 arrowline(c(xx[2]-ww[2]/2,xx[1]+ww[1]/2),c(yy,yy),col=amber)
 arrowline(c(xx[2]+ww[2]/2,xx[3]-ww[3]/2),c(yy,yy),col=amber)
 arrowline(c(xx[4]-ww[4]/2,xx[3]+ww[3]/2),c(yy,yy),col=amber)
 arrowline(c(xx[4]+ww[4]/2,xx[5]-ww[5]/2),c(yy,yy),col=amber)
 labels<-c('FT4','M','Creatinine','True GFR','CysC')
 for(i in seq_along(xx))box(xx[i],yy,ww[i],12,labels[i],unmeasured=i%in%c(2,4),adjusted=i==3)
 txt('Conditioning on creatinine can open this path if M remains unblocked.',85,151,7.5)
 txt('Dashed boxes: not fully measured. Dashed arrow: candidate pathway.',85,159,7.5)
 txt('Arrows express assumptions, not estimated effects or the cause of the sign reversal.',85,166,7.2)
 popViewport()
}
p<-file.path(out,'figures/FigureS3_hypothesized_DAG');w<-170/25.4;h<-170/25.4
svglite::svglite(paste0(p,'.svg'),width=w,height=h);draw();dev.off()
cairo_pdf(paste0(p,'.pdf'),width=w,height=h,family=font);draw();dev.off()
ragg::agg_tiff(paste0(p,'_600dpi.tiff'),width=w,height=h,units='in',res=600,compression='lzw',background='white');draw();dev.off()
ragg::agg_png(paste0(p,'_preview.png'),width=w,height=h,units='in',res=180,background='white');draw();dev.off()
pdftools::pdf_convert(paste0(p,'.pdf'),format='png',dpi=180,filenames=file.path(out,'qa/FigureS3_PDF.png'),verbose=FALSE)
fonts<-pdftools::pdf_fonts(paste0(p,'.pdf'));stopifnot(all(fonts$embedded))
pdftext<-pdftools::pdf_text(paste0(p,'.pdf'));stopifnot(grepl('True GFR',pdftext),grepl('noncausal backdoor',pdftext))
write.csv(edges,file.path(out,'evidence/DAG_edges.csv'),row.names=FALSE)
write_json(list(width_mm=170,height_mm=170,font_pt=c(7.2,9),all_fonts_embedded=all(fonts$embedded),acyclic=TRUE,conditioned_node='Cr',highlighted_path='FT4 <- M -> Cr <- GFR -> CysC',effect_sizes_estimated=FALSE,exploratory_post_result=TRUE,export_dpi=600),file.path(out,'evidence/DAG_QA.json'),auto_unbox=TRUE,pretty=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,'evidence/R_DAG_sessionInfo.txt'))

}

draw_associations(root)
draw_sensitivity(root)
draw_dag(root)
