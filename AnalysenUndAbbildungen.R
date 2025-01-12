###### Laden bzw Installieren der nötigen Pakete

################################################################################
# Liste der benötigten Pakete
required_packages <- unique(c(
  "tidyverse", "dplyr", "ggplot2", "rjson", "jsonlite", 
  "leaflet", "RCurl", "sf", "rgdal", "INLA", "INLAutils", 
  "cowplot", "pals", "viridisLite", "colorspace", "gganimate", 
  "gstat", "ggfortify", "forecast", "stars", "openxlsx", 
  "epm", "raster", "ggregplot", "nngeo"
))

# Funktion, um fehlende Pakete zu installieren und alle zu laden
install_and_load_packages <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      message(paste("Paket", pkg, "ist nicht installiert. Installation läuft..."))
      # Besondere Installation für INLA
      if (pkg == "INLA") {
        install.packages("INLA", repos = "https://inla.r-inla-download.org/R/stable")
      } else {
        install.packages(pkg, dependencies = TRUE)
      }
    }
    library(pkg, character.only = TRUE)
    message(paste("Paket", pkg, "wurde geladen."))
  }
}

# Aufruf der Funktion
install_and_load_packages(required_packages)

###############################################################################
###############################################################################
###### Vorbereitungen für die Analyse
###############################################################################

### Einladen der Funktion für die Modellselsktion vom GitHub Repositoriy
source(url("https://raw.githubusercontent.com/EichenbergBEF/AmphibienSachsen/refs/heads/main/Funktion_Modellselektion.R"))


## Url für Beispieldatensatz auf GitHub Repositoriy angeben
url <- "https://raw.githubusercontent.com/EichenbergBEF/AmphibienSachsen/main/Beispieldaten.csv"

MyData<- read.csv(url, stringsAsFactors = T)

#str(MyData)

###### Spezifizieren der Erwartungswerte (Priors) für die räumliche Komponente
## Grundlage: recherchen zu den Aktionsradien der Tiere
## im Beispieldatensatz wurden diese frei erfunden
## Prior: Erwartungswert
## alpha_scaled: erwartete Streuung um den Erwartungswert, skaliert zwischen 0.1 und 0.4
## nach der Spannweite der in der Literatur angegebenen Aktionsradien.
## kleinere Werte ind er Erwarteten Streuung bedeuten informiertere Priors

action_reprod<- data.frame("Art" = unique(MyData$Art),
                           "Species" = unique(MyData$Species),
                           "Spec_short" = unique(MyData$Spec_short),
                           "prior" = c(3500),
                           "alpha_scaled" = c(0.35)
                           )


################## räumliche Informationen für die Modelle vorbereiten

## Shapefiles für das Spatial Field von GitHub repository laden 

url2<- "https://raw.githubusercontent.com/EichenbergBEF/AmphibienSachsen/main/sachsenshape.rds"


sachsen<- readRDS(url(url2))

##### Definieren eines Working Directories
### Dieses müssen Sie selbst definieren. Hier werden alle weiteren Daten erzeugt bzw abgelegt

setwd(choose.dir())

################################################################################
################# Durchführen der Analyse mit dem Modelldatensatz ##############
################################################################################


for (i in unique(MyData$Spec_short)){
  
  print(paste0("Running Model selection for ",i," and calculating best Model"))
  print("############################################")
  print("WARNING: This can take very long for certain models.")
  print("Be assured, that the procedure works, even if you do not get response for quite some time.")
  print("In case of Problem,s an ERROR is thrown and the computation stops")
  print("############################################")

  # Erzeugen eines Ordners für den Output der Modelle
  dir.create(paste0(getwd(),"/Models"))

  
  dat<- MyData[MyData$Spec_short==i,]
  #dat <- dat %>% drop_na(Tstab_Summer_Vorj)
  #### Mesh machen
  
  Locs_spec<- cbind(unique(dat$X),unique(dat$Y))
  
  
  range_guess<- round(as.numeric(quantile(dist(Locs_spec, upper=T, diag=F)/1000,0.15)),0)
  MaxEdge<- range_guess/5
  
  mesh_spec<- inla.mesh.2d(boundary = sachsen,
                           loc = Locs_spec,
                           max.edge = c(1,5)*(MaxEdge*1000),
                           cutoff= (MaxEdge*1000)/5)
  
  ## SPDE machen
  
  #### HIer kann dann die Aktivitätsreichweite der entsprechenden Art rein
  pc.range<- action_reprod$prior[action_reprod$Spec_short == i]
  alpha.range<- action_reprod$alpha_scaled[action_reprod$Spec_short == i]
  
  spde_amphi<- inla.spde2.pcmatern(mesh_spec, 
                                   prior.range=c(pc.range*1000, alpha.range), 
                                   prior.sigma=c(0.5,0.5)) 
  
  A_amphi  <- inla.spde.make.A(mesh_spec, 
                               loc = cbind(dat$X,dat$Y))
  
  w.amphi<- inla.spde.make.index(name="w",
                                 n.spde = spde_amphi$n.spde)
  
  ####### den Stack machen
  
  N<- nrow(dat)
  modeldat<- data.frame("Intercept" = rep(1,N),
                        "Anzahl" = dat$Anzahl,
                        "Jahr" = dat$Jahr,
                        "Indiv" = dat$Anzahl,
                        "Jahr.1" = dat$Jahr-min(dat$Jahr)+1,
                        "Laenge" = MyStd(dat$Laenge),
                        "Dauer" = factor(dat$Dauer),
                        "Material" = factor(dat$Material),
                        "Temp_mean_Fruehling" = MyStd(dat$Temp_mean_Fruehling),
                        "Temp_stab_Fruehling" = MyStd(dat$Temp_stab_Fruehling),
                        "KWB_Fruehling" = MyStd(dat$KWB_Fruehling),
                        "Temp_mean_Fruehling_Vorj" = MyStd(dat$Temp_mean_Fruehling_Vorj),
                        "Temp_mean_Herbst_Vorj" = MyStd(dat$Temp_mean_Herbst_Vorj),
                        "Temp_mean_Winter_Vorj" = MyStd(dat$Temp_mean_Winter_Vorj),
                        "KWB_Fruehling_Vorj" = MyStd(dat$KWB_Fruehling_Vorj),
                        "KWB_Sommer_Vorj" = MyStd(dat$KWB_Sommer_Vorj),
                        "KWB_Herbst_Vorj" = MyStd(dat$KWB_Herbst_Vorj),
                        "KWB_Winter_Vorj" = MyStd(dat$KWB_Winter_Vorj)
                        )
  
  
  
  resp<- c("Anzahl")
  preds<- c("Laenge" , "Dauer", "Material",  
            "Temp_mean_Fruehling", "Temp_stab_Fruehling",
            "KWB_Fruehling",
            "Temp_mean_Fruehling_Vorj",
            "Temp_mean_Herbst_Vorj",
            "Temp_mean_Winter_Vorj",
            "KWB_Fruehling_Vorj",
            "KWB_Sommer_Vorj",
            "KWB_Herbst_Vorj",
            "KWB_Winter_Vorj")
  
  
  randoms<- c("Jahr.1")
  randomMods<- c("rw1")
    
  ### Die Antwortvariable ist nbinomial verteilt; daher muss die Zuordnung in INLA berücksichtigt werden
  family_funct<- c("nbinomial")

 ### Datenstack für INLA machen 
  Stack_amphi<- inla.stack(tag="Fit_amphi",
                           data= list(y=modeldat$Anzahl),
                           A = list(A_amphi,1),
                           effects = list(w.amphi,
                                          modeldat)
  )
  
  
  ## Modellselektion durchführen
  ## ACHTUNG: hierfür ist es notwendig, die Funktion einzulesen
  ## s. oben
  HostModelSel<- INLAModelSel_amphi(Response = resp,
                             Explanatory = preds,
                             Random = randoms,
                             RandomModel = randomMods,
                             Data = modeldat, 
                             Family= family_funct,
                             ScaleVariables = F,
                             Beep=F
                             )
 
  ### Extrahieren der finalen Covariablen für das optimale Modell
  Finalcovar <- HostModelSel$Removed[[length(HostModelSel$Removed)]]
  
  ### Formulieren der Modellformel
  ## hier ist wieder darauf zu achten, dass die räumliche
  ## und zeitliche Komponente explizit mit aufgenommen werden
  f_model<- as.formula(paste0("y", " ~ -1 + Intercept +",
                              paste(Finalcovar [!(Finalcovar %in% c("-1","Intercept"))], collapse = " + "),
                              "+ f(Jahr.1, model = 'rw1') + f(w, model = spde_amphi)"))
  
  ## Laufen lassen des optimalen Modells
  ## zur post-hoc bestimmung der Modellgüte
  ## muss hier control.compute(config = T) gesetzt werden!
  model_amphi<- inla(f_model,
                     family = family_funct,
                     data = inla.stack.data(Stack_amphi),
                     control.predictor = list(A = inla.stack.A((Stack_amphi)), compute = T, quantiles = c(0.025,0.975),link = 1),
                     control.compute = list(dic = TRUE, config =T, return.marginals.predictor=TRUE),
                     verbose = T
                     )
  
  
  saveRDS(file = paste0(getwd(),"/Models/best_model_",i,"_",family_funct,".rds"),object = model_amphi)
  
}


#### Durchführen der Modellgüte Analysen
### Goodness-of-fits werden anhand der Kuongruenz zwischen den beobachteten 
### Daten und den aufgrund der im Optimalen Modell identifizieren Paramer 
### Simulierten Daten ermittelt

### Als Metriken dienen
## a) Hellinger Distanz
## b) Kullbak-Leibler Divergenz

#### Funktionen defineiren
## Hellinger Distances
hellinger_distance <- function(obs, pred) {
  # Ensure inputs are valid
  if (length(obs) != length(pred)) {
    stop("The observed and predicted vectors must have the same length.")
  }
  
  # Normalize observed and predicted vectors to probabilities
  p <- obs / sum(obs)
  q <- pred / sum(pred)
  
  # Compute Hellinger distance
  distance <- sqrt(sum((sqrt(p) - sqrt(q))^2) / 2)
  
  return(distance)
}


## Kullback Leibler divergence
kl_divergence <- function(obs, pred) {
  # Ensure inputs are valid
  if (length(obs) != length(pred)) {
    stop("The observed and predicted vectors must have the same length.")
  }
  
  # Normalize observed and predicted vectors to probabilities
  p <- obs / sum(obs)
  q <- pred / sum(pred)
  
  # Replace zero probabilities in q with a small value to avoid division by zero
  q[q == 0] <- 1e-10
  
  # Compute KL divergence
  divergence <- sum(p * log(p / q), na.rm = TRUE)
  
  return(divergence)
}



#### Für Bayesian goodness-of-fit

### Tabelle vorbereiten, die die Werte für die Modellgüte
### aufnehmen

speclist<- unique(MyData$Spec_short)

GoodnessOfFits<- data.frame("Spec_short" = unique(MyData$Spec_short),
                            "Species" =unique(MyData$Species),
                            "Art" = unique(MyData$Art),
                            "HellingersDistance" = NA,
                            "KullbackLeibler" = NA)



for(i in speclist){
  
  ### Vorbereitung für Loop
  Species <- i
  Taxon <- GoodnessOfFits$Species[GoodnessOfFits$Spec_short == i]
  
  ### Erzeugen eines Unterordners zum Ablegen der Informationen
  dir.create(file.path(paste0(getwd(),"/Model_output/",i)),recursive=T)
  
  
  ### Optimales Modell laden
  mod<- readRDS(paste0(getwd(),"/Models/best_model_",i,"_nbinomial.rds"))
  
  
  GoodnessOfFits$Range_spat[which(GoodnessOfFits$Spec_short == i)]<- mod$summary.hyperpar$mean[grep("Range",rownames(mod$summary.hyperpar))]/1000
  GoodnessOfFits$sd_Spat[which(GoodnessOfFits$Spec_short == i)]<- mod$summary.hyperpar$sd[grep("Range",rownames(mod$summary.hyperpar))]/1000
  
  
  GoodnessOfFits$PC.Prior.U <- action_reprod$prior[match(GoodnessOfFits$Spec_short,action_reprod$Spec_short)]
  GoodnessOfFits$PC.Prior.prec <- action_reprod$alpha_scaled[match(GoodnessOfFits$Spec_short,action_reprod$Spec_short)]
  
  
  print(paste0("running goodness-of-fit assessment for ",i))
  
  ### Ouput vorbereiten
  filename<-paste0(getwd(),"/Model_output/",i)
  
  
  
  dat<- MyData[MyData$Spec_short==i,c("Anzahl","X","Y")]
  #### Mesh machen
  
  Locs_spec<- cbind(unique(dat$X),unique(dat$Y))
  
  
  range_guess<- round(as.numeric(quantile(dist(Locs_spec, upper=T, diag=F)/1000,0.15)),0)
  MaxEdge<- range_guess/5
  
  mesh_spec<- inla.mesh.2d(boundary = sachsen,
                           loc = Locs_spec,
                           max.edge = c(1,5)*(MaxEdge*1000),
                           cutoff= (MaxEdge*1000)/5)
  
  ## SPDE machen
  
  A_amphi  <- inla.spde.make.A(mesh_spec, 
                               loc = cbind(dat$X,dat$Y))
  
  
  y <- dat$Anzahl
  
  ##############################################################################
  ## Ziehen der Posterior Predictive Samples
  ##############################################################################
  
  n_samples <- 1000
  posterior_samples <- inla.posterior.sample(n_samples, mod)
  
  #### Auslesen der Modellparameter des optimnalen Modells
  modelparams<- unlist(strsplit(as.character(mod$.args$formula), "\\s+"))
  modelparams<- modelparams[!(modelparams %in% c("+","~","y","-1","+f(Jahr.1,","f(Jahr.1,",
                                                 "model","=","\"rw1\")","f(w,",
                                                 "spde_amphi)"))]
  
  MyParams<- c(modelparams,"Jahr.1")
  
  ####### Daten aus dem Optimalen modell einlesen, um Model fit zu ermöglichen
  simdata<- data.frame(mod$.args$data[which(names(mod$.args$data) %in% MyParams)])
  
  ## Formatieren der Daten, um sicherzugehen
  simdata$Jahr.1<- factor(simdata$Jahr.1)
  
  ## Model matrix für feste Variablen definieren
  modeldata<- as.matrix(model.matrix(as.formula(paste0("~ ",paste0(modelparams, collapse=" + "), " - 1")), data=simdata))
  
  ## Model Matrix für RW1 trend
  modelmatrix_temp<- model.matrix(~Jahr.1 -1, data=simdata)
  
  
  
  ##### Mittelwerte der linearen Prädiktioren über die Anzahl an Samples erzeugen 
  response_means <- lapply(1:length(posterior_samples), function(l) {
    sample <- posterior_samples[[l]]  # Access the i-th posterior sample
    
    ## Extrahiere Betas für Fixed Effects
    rownums_fixed<- unlist(as.vector(lapply(modelparams, function(x)
      grep(x,rownames(sample$latent), fixed =T))))
    rownums_fixed<- rownums_fixed[!duplicated(rownums_fixed)]
    Betas_fixed <- sample$latent[rownums_fixed]
    
    ### Berechne etas für FixedPart
    FixedPart<- modeldata %*% Betas_fixed
    
    
    ### Extrahiere Betas für RW1 Trend
    rownums_temp<- as.vector(sapply("Jahr", function(x)
      grep(x,rownames(sample$latent), fixed =T)))
    
    Betas_temp <- c(sample$latent[rownums_temp])
    
    
    ## Berechne etas für RW1 Trend
    TempPart<- modelmatrix_temp %*%  Betas_temp
    
    
    ### Extrahiere etas für röumlichen Part
    rownums_spat<- as.vector(sapply("w", function(x)
      grep(x,rownames(sample$latent), fixed =T)))
    
    Betas_spat<- c(sample$latent[rownums_spat])
    ### Berechne räuml. Effekt
    
    SpatPart<- A_amphi %*%  Betas_spat
    
    mu<- round(exp(FixedPart + TempPart + SpatPart),0)
    mu<- as.vector(mu)
    return(mu)
    
  })
  
  
  # Extrahieren der simulierten Overdispersion Parameter (phi)
  phi_samples <- sapply(1:length(posterior_samples), function(k) {
    sample <- posterior_samples[[k]]
    phi<- sample$hyperpar[grep("nbinomial", names(sample$hyperpar))]
    #phi <- exp(log_phi)
    return(phi)
  })
  
  
  
  
  # Simulieren der Vorhergesagten Anzahl
  simulated_responses <- lapply(1:n_samples, function(i) {
    rnbinom(length(response_means[[i]]), size = phi_samples[i], mu = response_means[[i]])
  })
  
  
  ### Extraktion der observed data (ohne NA's)
  observed_mean<- y
  valid_indices <- which(!is.na(y))  # Indices of non-NA entries in y
  y_clean <- y[valid_indices]        # Filtered y (observed values)
  y_pred_clean <- round(observed_mean[valid_indices],0)
  
  ### Berechnen der goodnesses-of-fits:
  
  ### Überlapp der Oberved vs mean posterior predictive data
  ### Mittelwert der Positerior distributions erzeugen
  sim_pred_data<- colMeans(do.call(rbind, simulated_responses))[valid_indices]
  
  
  
  #### Hellinger Distance
  
  hell_dist<-hellinger_distance(y_clean,sim_pred_data)
  kld<- kl_divergence(y_clean,sim_pred_data)
  
  GoodnessOfFits$HellingersDistance[which(GoodnessOfFits$Spec_short == i)] <- hell_dist
  GoodnessOfFits$KullbackLeibler [which(GoodnessOfFits$Spec_short == i)] <- kld
  
  ## Plot der Überlappung machen
  plot_data <- data.frame(
    value = c(y_clean, sim_pred_data),
    group = c(rep("Beobachtet", length(y_clean)), 
              rep("Vorhergesagt", length(sim_pred_data)))
  )
  
  gof_plot<- ggplot(plot_data, aes(x = value, color = group, fill = group)) +
    geom_density(alpha = 0.5) +
    scale_color_manual(values = c("darkgray", "salmon")) +  # Line colors
    scale_fill_manual(values = c("cyan", "salmon")) +   # Fill colors
    labs(title = paste0("Beobachtete Anzahl vs simulierte Anzahl für die Art ", Taxon),
         subtitle = paste0("Gemittelt über ", n_samples, " Simulationen"),
         x = "Anzahl",
         y = "Wahrscheinlichkeitsdichte",
         color = "Data Type",
         fill = "Data Type") +  # Legend title for both color and fill
    theme(plot.title = element_text(hjust = 0.5),  # Center the title
          plot.subtitle = element_text(hjust = 0.5)  # Center the subtitle
    ) +
    annotate(
      "text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.05, 
      label = paste0("Hellinger Dist: ", round(hell_dist, 3), 
                     "\nKL Divergence: ", round(kld, 3)),
      size = 4, color = "black"
    )
  
  ## Abspeichern des grafischen Outputs
  ggsave(filename = paste0(getwd(),"/Model_output/",i,"/GoodnessOfFit_",i,".png"), device = "png", width = 1754, height= 1240, dpi = 150, units = "px")
  
  
  ## Abspeichern des nummerischen Ouitputs
  write.xlsx(GoodnessOfFits, file= paste0(getwd(),"/Model_output/GoodnessOfFit_table.xlsx"))
  
}



################################################################################
######## Abbildungen aus den Modellparametern eruegen ##########################
################################################################################

for(i in unique(MyData$Spec_short)){
  mod<- readRDS(paste0(getwd(),"/Models/best_model_",i,"_nbinomial.rds"))
  
  
  Art<- unique(na.omit(MyData$Art[which(MyData$Spec_short==i)]))
  Species<- unique(na.omit(MyData$Species[which(MyData$Spec_short==i)]))
  
  filename<-paste0(getwd(),"/Model_output/",i)
  
  
  
  sachsen_utm<- spTransform(sachsen, CRS("+proj=utm +zone=33 +ellps=WGS84 +units=m +no_defs"))
  sachsen_sf<- st_as_sf(sachsen_utm)
  
  
  #### Summary-table dr Modellparameter machen
  model_summary<- as.data.frame(summary(mod)$fixed)
  
  model_summary$meaningful<-c("")
  model_summary$meaningful[which(sign(model_summary$`0.025quant`)==sign(model_summary$`0.975quant`))]<- c("Ja")
  
  names(model_summary) <- c("Mittelwert","sd","2.5% Quant", "50% Quant", "97.5% Quant", "Modus", "kld","Bedeutsam")
  model_summary$Parameter <- row.names(model_summary)
  write.xlsx(model_summary[,c(9,1,2,3,5,8)], file = paste0(filename,"/summary_model_fixed.xlsx"))
  
  
  ## Vorbereiten Datensatz für Abbildungen
  Spec_dat<- MyData[which(MyData$Spec_short == i),]
  
  
  Spec_dat$preds<- mod$summary.fitted.values$mean[1:nrow(Spec_dat)]
  Spec_dat$upper<- mod$summary.fitted.values$`0.975quant`[1:nrow(Spec_dat)]
  Spec_dat$lower<- mod$summary.fitted.values$`0.025quant`[1:nrow(Spec_dat)]
  
  
  ##### Abbldung für Vorhersagen an den jeweiligen Zaunstandorten
  
  predictions_at_locs<- ggplot() + geom_point(data= Spec_dat, aes(x=Jahr,y=Anzahl), col="gray50",alpha=0.6) + 
    geom_point(data= Spec_dat, aes(x=Jahr,y=preds),col="red")+
    geom_line(data= Spec_dat, aes(x=Jahr,y=preds),col="black", lty=2) + 
    geom_ribbon(data=Spec_dat, aes(x = Jahr, 
                                   ymax = upper, 
                                   ymin = lower),
                fill = grey(0.5),
                alpha = 0.4) + 
    facet_wrap(~Standort,scales="free_y") + theme(plot.title = element_text(hjust=0.5,face = "italic")) +
    ggtitle(paste(Art,Species,sep=" / "))  
  
  ggsave(predictions_at_locs,file = paste0(filename,"/Vorhersagen_Standorte.png"), device = "png", width = 1754, height= 1240, dpi = 150, units = "px")
  
  ## Beobachtet vs. Predicted
  
  beobachtet_predicted<- ggplot() + geom_point(data= Spec_dat, aes(x=Anzahl,y=preds)) + 
    geom_pointrange(data= Spec_dat, aes(x=Anzahl,y=preds,ymin = lower, ymax = upper),col="gray60")+
    xlab("Beobachtet") + ylab("Vorhersage") + geom_smooth(data=Spec_dat, aes(x=Anzahl, y=preds),method="lm") + 
    theme(plot.title = element_text(hjust=0.5,face = "italic")) +
    ggtitle(paste(Art,Species,sep=" / "))
  
  ggsave(beobachtet_predicted,file = paste0(filename,"/Observed_Fitted.png"), device = "png", width = 1754, height= 1240, dpi = 150, units = "px")
  
  
  ##### Zeitlicher Trend
  trend_dat<- data.frame(Jahr=sort(unique(Spec_dat$Jahr)),
                         mean= mod$summary.random$Jahr.1$mean,
                         upper=mod$summary.random$Jahr.1$`0.975quant`,
                         lower=mod$summary.random$Jahr.1$`0.025quant`)
  
  
  trend<-  ggplot() + geom_point(data=trend_dat,aes(x=Jahr,y=mean),col="red") +
    geom_line(data=trend_dat,aes(x=Jahr,y=mean),col="black") + 
    geom_hline(yintercept = 0, lty=2)+
    geom_ribbon(data=trend_dat, aes(x = Jahr, 
                                    ymax = upper, 
                                    ymin = lower),
                fill = grey(0.5),
                alpha = 0.4)+ theme(plot.title = element_text(hjust=0.5,face = "italic"))+
    theme(axis.title.y = element_text(size = 14),
          axis.title.x = element_text(size = 14)) +
    ylab("Einfluss der \n Vorjahrespopulationsgröße") +
    ggtitle(paste(Art,Species,sep=" / "))    
  
  ggsave(trend,file = paste0(filename,"/timetrend.png"), device = "png", width = 1754, height= 1240, dpi = 150, units = "px")
  
  
  
  ######### Für Karten des räumlichen Effektes
  
  ##### Für Karten mit Mittelwert und Standardabweichung vom Mittleren Trend
  
  locs<- cbind(Spec_dat$X,Spec_dat$Y)
  
  loc_range<- cbind(Spec_dat$X[!duplicated(Spec_dat$Standort)],Spec_dat$Y[!duplicated(Spec_dat$Standort)])
  exp.range<- round(as.numeric(quantile(dist(loc_range, upper=T, diag=F)/1000,0.15)),0)
  MaxEdge <- exp.range/5
  
  ConvHull<- inla.nonconvex.hull(locs)
  
  mesh<-inla.mesh.2d(boundary = sachsen,loc=locs,
                     max.edge=c(1,5)*(MaxEdge*1000),
                     cutoff= (MaxEdge*1000)/5) 
  
  
  x.grid <- seq(min(st_coordinates(sachsen_sf)[,1]),max(st_coordinates(sachsen_sf)[,1]), length.out=300)
  y.grid <- seq(min(st_coordinates(sachsen_sf)[,2]),max(st_coordinates(sachsen_sf)[,2]), length.out=300)
  
  Grid <- expand.grid(X = x.grid, 
                      Y = y.grid)
  
  
  wproj<-  inla.mesh.projector(mesh, loc = as.matrix(Grid))
  w.pm.m<- inla.mesh.project(wproj,mod$summary.random$w$mean)
  w.pm.sd<- inla.mesh.project(wproj,mod$summary.random$w$sd)
  
  
  Grid$w.pm <- as.vector(w.pm.m)     
  Grid$w.sd <- as.vector(w.pm.sd)               
  
  GridNoNA  <- na.exclude(Grid)
  
  GridNoNA<- st_as_sf(GridNoNA, coords=c("X","Y"), crs=st_crs(sachsen_sf))
  
  GridNoNA<- st_intersection(GridNoNA,sachsen_sf)
  
  GridNoNA_ras_mean <- st_rasterize(GridNoNA %>% dplyr::select(w.pm, geometry))
  GridNoNA_ras_mean<- st_crop(GridNoNA_ras_mean,sachsen_sf)
  
  GridNoNA_ras_sd<- st_rasterize(GridNoNA %>% dplyr::select(w.sd, geometry))
  GridNoNA_ras_sd<- st_crop(GridNoNA_ras_sd,sachsen_sf)
  
  saveRDS(GridNoNA_ras_mean,paste0(filename,"/GridNoNA_ras_mean_",i,".rds"))
  saveRDS(GridNoNA_ras_sd,paste0(filename,"/GridNoNA_ras_sd_",i,".rds"))
  
  
  plot_spat_mean<- ggplot() + geom_sf(data= GridNoNA, aes(col=w.pm),shape=15, size=2)+ geom_point(data=Spec_dat, aes(x=X, y=Y), col="gray60", alpha=0.04) +
    scale_color_continuous_divergingx(palette = "RdBu", mid=0,name = "Mittelw. des räumlichen \n Einflusses") + 
    xlab("Rechtswert") + ylab("Hochwert") + 
    theme(plot.title = element_text(hjust=0.5), 
          axis.text.x=element_text(angle=23, hjust=1,size = 9),
          axis.text.y=element_text(angle=23, hjust=1,size = 9),
          axis.title = element_blank(),
          strip.text.x = element_text(size = 12,face="bold.italic"),
          strip.text.y = element_text(size = 12,face="bold.italic"),
          legend.title=element_text(hjust=0.5, vjust=0.5),
          legend.text = element_text(angle=30, hjust=1),
          legend.position="right")+
    ggtitle(paste(Art,Species,sep=" / ")) +   
    geom_sf(data = sachsen_sf, fill = NA, colour = "black", size = 0.8)+ theme(plot.title = element_text(hjust=0.5,face = "italic"))
  
  ggsave(plot_spat_mean,file = paste0(filename,"/Spat_Resid_mean.png"), device = "png", width = 1754, height= 1240, dpi = 150, units = "px")
  
  plot_spat_sd<- ggplot() + geom_sf(data= GridNoNA, aes(col=w.sd),shape=15, size=2)+ geom_point(data=Spec_dat, aes(x=X, y=Y), col="gray60", alpha=0.04) +
    scale_color_continuous_sequential(palette = "Blues 3",name = "SD des räumlichen \n Einflusses") + 
    xlab("Rechtswert") + ylab("Hochwert") + 
    theme(plot.title = element_text(hjust=0.5), 
          axis.text.x=element_text(angle=23, hjust=1,size = 9),
          axis.text.y=element_text(angle=23, hjust=1,size = 9),
          axis.title = element_blank(),
          strip.text.x = element_text(size = 12,face="bold.italic"),
          strip.text.y = element_text(size = 12,face="bold.italic"),
          legend.title=element_text(hjust=0.5, vjust=0.5),
          legend.text = element_text(angle=30, hjust=1),
          legend.position="right")+
    ggtitle(paste(Art,Species,sep=" / ")) +
    geom_sf(data = sachsen_sf, fill = NA, colour = "black", size = 0.8)+ theme(plot.title = element_text(hjust=0.5,face = "italic"))
  
  ggsave(plot_spat_sd,file = paste0(filename,"/Spat_Resid_sd.png"), device = "png", width = 1754, height= 1240, dpi = 150, units = "px")
  
}
