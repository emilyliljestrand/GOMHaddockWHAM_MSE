# Gulf of Maine Haddock MSE in WHAM
# Emily Liljestrand
# Created: Oct 25, 2023

# Remove Old Objects and set working directory to file location:
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
originalwd = getwd()

#Install and load packages:
get_F_from_catch <- function(om, year, catch, Finit = 0.1, maxF = 10){
  rep = om$report()
  #print(year)
  naa = rep$NAA[year,]
  Maa = rep$MAA[year,]
  sel_tot = rep$FAA_tot[year,]/max(rep$FAA_tot[year,])
  waa = om$input$data$waa[om$input$data$waa_pointer_totcatch, year,]
  get_catch = function(log_F, naa, sel, waa, Maa){
    Faa = exp(log_F) * sel_tot
    Zaa = Maa + Faa
    Catch = 0
    for(a  in 1:length(naa)) Catch = Catch + waa[a] * naa[a] * Faa[a] *(1 - exp(-Zaa[a]))/Zaa[a];
    return(Catch)
  }
  obj = function(log_F) (catch - get_catch(log_F, naa, sel_tot, waa, Maa))^2
  opt = try(nlminb(log(Finit), obj))
  if(!is.character(opt)) Fsolve = exp(opt$par)[1] else Fsolve = maxF
  if(Fsolve>10) Fsolve = maxF
  print(paste0("Fsolve: ", Fsolve))
  return(Fsolve)
}

# install.packages('TMB', type = 'source')
# devtools::install_github("emilylil/wham")
# devtools::install_github("timjmiller/wham", dependencies=TRUE)
# devtools::install_github("timjmiller/wham", dependencies=TRUE,ref='dvel')
# devtools::install_github("timjmiller/wham", dependencies=TRUE, ref = "master", lib = "C:/Users/emily.liljestrand/AppData/Local/Programs/R/R-4.3.1/library/wham", INSTALL_opts=c("--no-multiarch"))
# devtools::install_github("timjmiller/wham", dependencies=TRUE, ref = "lab", lib = "C:/Users/emily.liljestrand/AppData/Local/Programs/R/R-4.3.1/library/multi_wham", INSTALL_opts=c("--no-multiarch"))
# install.packages(c("wham","TMB"))

library(wham, lib.loc = "C:/Users/emily.liljestrand/AppData/Local/Programs/R/R-4.3.1/library/wham")
library(TMB)
library(dplyr)
library(ggplot2)

##################################### DATA #####################################

#Read in Gulf of Maine Haddock Data from ASAP Dat File
# GOM_HADDOCK_DAT <- read_asap3_dat("GOM_HADDOCK_ASAP_2021_BASE_NEWCALIB.DAT")
#Read in Gulf of Maine Haddock Data from ASAP Dat File, version for selectivity random effects (one block for all catch)
GOM_HADDOCK_DAT <- read_asap3_dat("GOM_HADDOCK_ASAP_2021_1BLOCK_NEWCALIB.DAT")

#Read in existing seeds file, or make one if it doesn't exist
if(!file.exists("SIM.seeds.csv"))
{
  set.seed(799291)
  nseeds <- 1000
  r.seed.set <- trunc(1e7*runif(n=nseeds),7) + trunc(1e3*runif(n=nseeds),3)
  write.table(r.seed.set, file="SIM.seeds.csv", quote=F, row.names=F, col.names=F, sep=",")
} else(r.seed.set <- read.table("SIM.seeds.csv"))

# ##################################### MODEL(s) #####################################
#
# ------------------------------------ OM with selectivity random effects ------------------------------------

# Set and fit the OM, if the file doesn't exist
input <- prepare_wham_input(GOM_HADDOCK_DAT, recruit_model=2, model_name="GOMHaddock_RW_Recruitment",
                            selectivity=list(model=rep("age-specific",3),
                                             re=c(rep('ar1_y',1),rep("none",2)),
                                             initial_pars <- list(c(0.005,0.08,0.25,0.5,0.6,0.75,1,1,1),
                                                                  c(0.2,0.4,0.8,1,1,1,1,1,1),
                                                                  c(0.1,0.3,0.5,0.8,0.9,1,1,1,1)),
                                             # fix_pars=c(vector(mode = "list", length = 1),list(c(4:9),c(6:9)))),
                                             fix_pars=c(list(c(7:9),c(4:9),c(5:9)))),
                            NAA_re=list(sigma='rec+1',cor='ar1_y'),
                            age_comp='multinomial')

# Mapping of Fishing Mortality
# input$par$log_F1 <- log(3)
# input$map$log_F1 <- as.factor(NA)
# input$map$F_devs <- as.factor(matrix(data=NA,nrow=42,ncol=1))

# Mapping the NAA in the first year
# input$map$log_N1_pars=as.factor(matrix(data=NA,nrow=1,ncol=9)) # Fix all values
# input$map$log_N1_pars=as.factor(matrix(data=c(1,2,3,4,5,6,NA,NA,7),nrow=1,ncol=9)) # Fix ages 7-8
input$map$log_N1_pars=as.factor(matrix(data=c(1,2,3,4,5,6,NA,NA,NA),nrow=1,ncol=9)) # Fix ages 7-8
input$par$log_N1_pars[7:8] <- c(log(5),log(1))
input$par$log_N1_pars[9] <- log(20)

# Fit model, if you haven't already fit this model:
if(!file.exists("mod.OM.rds"))
{
  mod.OM <- try(fit_wham(input, do.osa = T, do.check=T,do.retro=F)) 
  saveRDS(mod.OM, file='mod.OM.rds')
  mod.OM.proj <- project_wham(mod.OM,proj.opts=list(use.FXSPR=T))
  saveRDS(mod.OM.proj, file='mod.OM.proj.rds')

  Foldername <- "mod.OM.Output"
  if(!dir.exists(Foldername)) dir.create(Foldername)
  setwd(file.path(getwd(),Foldername))
  # plot_wham_output(mod=mod.OM,out.type="html")
  plot_wham_output(mod=mod.OM)
  setwd(originalwd)

} else {mod.OM<-readRDS('mod.OM.rds')
        mod.OM.proj <-readRDS('mod.OM.proj.rds')}
simpar_base <- mod.OM$env$last.par.best

# Beginning of MSE -------------------------------------------

# Vector of all the F reference points (e.g., F40%)
Fref <- c()
Fref <- c(Fref,0.3)
# Fref <- c(Fref,exp(last(mod.OM$rep$log_FXSPR)))

n_MSE_years <- 8

for(i in 1:n_MSE_years)
{
# Fit a "pseudo" model to have a WHAM model with the correct number of containers
# Specify data for the "pseudo" model with baseline data plus however many extra years of "pseudo" data
  GOM_HADDOCK_DAT_MSE <- GOM_HADDOCK_DAT
  GOM_HADDOCK_DAT_MSE$dat$n_years <-GOM_HADDOCK_DAT$dat$n_years + i
  GOM_HADDOCK_DAT_MSE$dat$M <- rbind(GOM_HADDOCK_DAT$dat$M,matrix(rep(mean(GOM_HADDOCK_DAT$dat$M),(i * GOM_HADDOCK_DAT$dat$n_ages)),nrow=i))
  GOM_HADDOCK_DAT_MSE$dat$maturity <- rbind(GOM_HADDOCK_DAT$dat$maturity,matrix(rep(apply(GOM_HADDOCK_DAT$dat$maturity,2,mean),i),nrow=i))
  GOM_HADDOCK_DAT_MSE$dat$WAA_mats[[1]] <- rbind(GOM_HADDOCK_DAT$dat$WAA_mats[[1]],matrix(rep(apply(GOM_HADDOCK_DAT$dat$WAA_mats[[1]],2,mean),i),nrow=i))
  GOM_HADDOCK_DAT_MSE$dat$WAA_mats[[2]] <- rbind(GOM_HADDOCK_DAT$dat$WAA_mats[[2]],matrix(rep(apply(GOM_HADDOCK_DAT$dat$WAA_mats[[2]],2,mean),i),nrow=i))
  GOM_HADDOCK_DAT_MSE$dat$sel_block_assign[[1]] <- c(GOM_HADDOCK_DAT$dat$sel_block_assign[[1]],rep(last(GOM_HADDOCK_DAT$dat$sel_block_assign[[1]]),i))
  GOM_HADDOCK_DAT_MSE$dat$CAA_mats[[1]] <- rbind(GOM_HADDOCK_DAT$dat$CAA_mats[[1]],matrix(rep(apply(GOM_HADDOCK_DAT$dat$CAA_mats[[1]],2,mean),i),nrow=i))
  GOM_HADDOCK_DAT_MSE$dat$DAA_mats[[1]] <- rbind(GOM_HADDOCK_DAT$dat$DAA_mats[[1]],matrix(rep(apply(GOM_HADDOCK_DAT$dat$DAA_mats[[1]],2,mean),i),nrow=i))
  GOM_HADDOCK_DAT_MSE$dat$prop_rel_mats[[1]] <- rbind(GOM_HADDOCK_DAT$dat$prop_rel_mats[[1]],matrix(rep(apply(GOM_HADDOCK_DAT$dat$prop_rel_mats[[1]],2,mean),i),nrow=i))
  GOM_HADDOCK_DAT_MSE$dat$IAA_mats[[1]] <- rbind(GOM_HADDOCK_DAT$dat$IAA_mats[[1]],matrix(rep(apply(GOM_HADDOCK_DAT$dat$IAA_mats[[1]],2,mean),i),nrow=i))
  GOM_HADDOCK_DAT_MSE$dat$IAA_mats[[2]] <- rbind(GOM_HADDOCK_DAT$dat$IAA_mats[[2]],matrix(rep(apply(GOM_HADDOCK_DAT$dat$IAA_mats[[2]],2,mean),i),nrow=i))
  GOM_HADDOCK_DAT_MSE$dat$recruit_cv <- c(GOM_HADDOCK_DAT$dat$recruit_cv,rep(last(GOM_HADDOCK_DAT$dat$recruit_cv),i))
  GOM_HADDOCK_DAT_MSE$dat$catch_cv <- rbind(GOM_HADDOCK_DAT$dat$catch_cv,matrix(rep(last(GOM_HADDOCK_DAT$dat$catch_cv),i),nrow=i))
  GOM_HADDOCK_DAT_MSE$dat$discard_cv <- rbind(GOM_HADDOCK_DAT$dat$discard_cv,matrix(rep(last(GOM_HADDOCK_DAT$dat$discard_cv),i),nrow=i))
  GOM_HADDOCK_DAT_MSE$dat$catch_Neff <- rbind(GOM_HADDOCK_DAT$dat$catch_Neff,matrix(rep(last(GOM_HADDOCK_DAT$dat$catch_Neff),i),nrow=i))
  GOM_HADDOCK_DAT_MSE$dat$discard_Neff <- rbind(GOM_HADDOCK_DAT$dat$discard_Neff,matrix(rep(last(GOM_HADDOCK_DAT$dat$discard_Neff),i),nrow=i))
  
  # Use the data and pseudo data to make a pseudo input with correct specification
  input_MSE <- prepare_wham_input(GOM_HADDOCK_DAT_MSE, recruit_model=2, model_name="GOMHaddock_RW_Recruitment",
                              selectivity=list(model=rep("age-specific",3),
                                               re=c(rep('ar1_y',1),rep("none",2)),
                                               initial_pars <- list(c(0.005,0.08,0.25,0.5,0.6,0.75,1,1,1),
                                                                    c(0.2,0.4,0.8,1,1,1,1,1,1),
                                                                    c(0.1,0.3,0.5,0.8,0.9,1,1,1,1)),
                                               # fix_pars=c(vector(mode = "list", length = 1),list(c(4:9),c(6:9)))),
                                               fix_pars=c(list(c(7:9),c(4:9),c(5:9)))),
                              NAA_re=list(sigma='rec+1',cor='ar1_y'),
                              age_comp='multinomial')
  input_MSE$map$log_N1_pars=as.factor(matrix(data=c(1,2,3,4,5,6,NA,NA,NA),nrow=1,ncol=9)) # Fix ages 7-8
  input_MSE$par$log_N1_pars[7:8] <- c(log(5),log(1))
  input_MSE$par$log_N1_pars[9] <- log(20)
  
  # Fit the pseudo model to have a basis for an OPERATING model
  mod.OM_MSE <- try(fit_wham(input_MSE, do.osa = F, do.check=F,do.retro=F))
  
  # Get the list of fitted parameters from pseudo model
  simpar <- mod.OM_MSE$env$last.par.best
  
  # Replace every value in the parameter object with the baseline values (from simpar_base)
  # Replace all the NEW parameters (from new years of data) with simujaltion values (F_devs, log_NAA, selpars_re)
  simpar[names(simpar)=="mean_rec_pars"] = simpar_base[names(simpar_base)=="mean_rec_pars"]
  simpar[names(simpar)=="logit_q"] = simpar_base[names(simpar_base)=="logit_q"]
  simpar[names(simpar)=="log_F1"] = simpar_base[names(simpar_base)=="log_F1"]
  simpar[names(simpar)=="F_devs"] = c(simpar_base[names(simpar_base)=="F_devs"],(log(Fref)-simpar_base[names(simpar_base)=="log_F1"])) #####Replace with F40% or whatever
  
  simpar[names(simpar)=="log_N1_pars"] = simpar_base[names(simpar_base)=="log_N1_pars"]
  simpar[names(simpar)=="log_NAA_sigma"] = simpar_base[names(simpar_base)=="log_NAA_sigma"]
  simpar[names(simpar)=="trans_NAA_rho"] = simpar_base[names(simpar_base)=="trans_NAA_rho"]
  
  simpar[names(simpar)=="log_NAA"] = c(simpar_base[names(simpar_base)=="log_NAA"],rep(rep(mean(simpar_base[names(simpar_base)=="log_NAA"]),GOM_HADDOCK_DAT$dat$n_ages),i))
  simpar[names(simpar)=="logit_selpars"] = simpar_base[names(simpar_base)=="logit_selpars"]
  simpar[names(simpar)=="selpars_re"] = c(simpar_base[names(simpar_base)=="selpars_re"],rep(0,i))
  simpar[names(simpar)=="sel_repars"] = simpar_base[names(simpar_base)=="sel_repars"]
  
  # Run Operating Model:
  simdata <- mod.OM_MSE$simulate(par=simpar,complete=T)
  
  # Set up the input for the estimation model by replacing catch/index with simulated values from OM:
  siminput <- input_MSE
  siminput$data$catch_paa <- simdata$catch_paa
  siminput$data$index_paa <- simdata$index_paa 
  siminput$data$agg_catch <- simdata$agg_catch
  siminput$data$agg_indices <- simdata$agg_indices
  siminput$data$obsvec <- simdata$obsvec
  siminput$data$index_Neff <-siminput$data$index_Neff 
  siminput$data$catch_Neff <-siminput$data$catch_Neff 
  siminput$data$agg_catch_sigma <- siminput$data$agg_catch_sigma
  siminput$data$agg_index_sigma <- siminput$data$agg_index_sigma
  
  # Fit estimation model
  simmod <- fit_wham(siminput,do.osa=F,do.check=F,do.retro=F)
  
  # Collect values of interest, get the recommended fishing mortality and close the loop
  # Fref <- c(Fref,exp(last(simmod$rep$log_FXSPR)))
  Fref <- c(Fref,0.3)
}

print(Fref)
print(simdata$SSB)