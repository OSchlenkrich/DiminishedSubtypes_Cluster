library(data.table)
library(fpc)
library(dplyr)
library(fclust)
library(parallel)
library(cluster)
library(ggplot2)
library(tidyr)
library(forcats)
library(gridExtra)
library(ggpubr)
library(scales)

# used by ggplot: N per cluster
stat_box_count <- function(y, upper_limit = 0.3) {
  return( 
    data.frame(
      y = 0.95 * upper_limit,
      label = paste('n =', length(y), '\n')
    )
  )
}


# Create Dataset ####
dmx_data = fread("Datasets/DemocracyMatrix_v2.csv", encoding="UTF-8") #  please replace the path

cluster_data = dmx_data %>%
  filter(classification_context == "Working Democracy" | classification_context == "Deficient Democracy") %>% 
  rename(freedom = freedom_dim_index_context,
         equality = equality_dim_index_context,
         control = control_dim_index_context) %>% 
  rowwise() %>%
  mutate(mean_dim =  mean(c(freedom, equality, control)),
         freedom = (freedom - mean_dim),
         equality = (equality - mean_dim),
         control = (control - mean_dim),
         
  ) %>%
  ungroup() %>% 
  select(mean_dim, freedom, equality, control) %>% 
  mutate_all(funs(uv = (. - min(., na.rm=T))/(max(., na.rm=T)-min(., na.rm=T)))) %>% 
  mutate(mean_dim_uv = mean_dim_uv * 0.5)  %>% 
  select_at(vars(ends_with("_uv"))) 



# Cluster Benchmark ####

interface_FKM = function(data, k, method) {
  FKM_results = FKM(X = data, k, maxit=2000, stand=0, RS=25)
  hc_classes = as.integer(FKM_results$clus[,1])
  #print(FKM_results$iter)
  
  cluster_sol = matrix(hc_classes, nrow=length(hc_classes), ncol=k)
  clusterlist =list()
  
  for (kk in 1:k) {
    clusterlist[[kk]] = if_else(cluster_sol[,kk]==kk, T, F)
  }
  
  make_list = list(
    result = hc_classes,
    nc=k,
    clusterlist = clusterlist, 
    partition=hc_classes,
    clustermethod="FKM"
  )
  
  return(make_list)
}

interface_FKMmed = function(data, k, method) {
  
  FKMmed_results = FKM.med(X = data, k, maxit=30, stand=0, RS=2)
  hc_classes = as.integer(FKMmed_results$clus[,1])
  print(FKMmed_results$iter)
  
  cluster_sol = matrix(hc_classes, nrow=length(hc_classes), ncol=k)
  
  clusterlist =list()
  
  for (kk in 1:k) {
    clusterlist[[kk]] = if_else(cluster_sol[,kk]==kk, T, F)
  }
  
  make_list = list(
    result = hc_classes,
    nc=k,
    clusterlist = clusterlist, 
    partition=hc_classes,
    clustermethod="FKMmed"
  )
  
  return(make_list)
}

# I modified the function, so that it calculates fewer cluster statistics. The function already inverses the statistics
# so that large values are good
source("Analyse/Cluster/clusterbenchstats_mod.R")
clustermethodpars <- list()
clustermethodpars[[4]] <- list()
clustermethodpars[[4]]$method <- "average"
clustermethodpars[[1]] <- list()
clustermethodpars[[1]]$method <- ""
clustermethodpars[[2]] <- list()
clustermethodpars[[2]]$runs <- 20
clustermethodpars[[3]] <- list()

distmethod <- c(F,F,F,F)

methodname <- c("FKM", "kmeans", "pam", "average")


set.seed(1234)
bootclassif <- c("centroid", "centroid", "centroid", "averagedist")
benchmark_results = clusterbenchstats_mod(cluster_data, G=2:10,
                                          diss = F,
                                          clustermethod =  c("interface_FKM","kmeansCBI", "pamkCBI", "hclustCBI"),
                                          clustermethodpars = clustermethodpars,
                                          distmethod = distmethod,
                                          methodnames = methodname,
                                          scaling=F,
                                          multicore = T,
                                          cores = 10,
                                          useboot = T,
                                          bootclassif=bootclassif,
                                          bootmethod="nselectboot",
                                          bootruns=100,
                                          nnruns = 100,
                                          fnruns = 100,
                                          avenruns = 100,
                                          kmruns = 100,
                                          useallg = F, 
                                          trace = F)





# Results Benchmark ####



aggregated_indizes = print(benchmark_results$sstat,
                           aggregate=TRUE,
                           weights=c(1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0))

bench_results_df = as.data.frame(do.call(rbind, benchmark_results$sstat[[1]])) %>% 
  mutate_all(funs(unlist)) %>% 
  mutate(method = "FCM",
         cluster = 1:10,
         nr = "Solution",
         aggregate = c(NA, unlist(aggregated_indizes[[17]][1,-1]))) %>% 
  bind_rows(as.data.frame(do.call(rbind, benchmark_results$sstat[[2]])) %>% 
              mutate_all(funs(unlist)) %>% 
              mutate(method = "kmeans",
                     cluster = 1:10,
                     nr = "Solution",
                     aggregate = c(NA, unlist(aggregated_indizes[[17]][2,-1])))) %>% 
  bind_rows(as.data.frame(do.call(rbind, benchmark_results$sstat[[3]])) %>% 
              mutate_all(funs(unlist)) %>% 
              mutate(method = "pam",
                     cluster = 1:10,
                     nr = "Solution",
                     aggregate = c(NA, unlist(aggregated_indizes[[17]][3,-1]))) ) %>% 
  bind_rows(as.data.frame(do.call(rbind, benchmark_results$sstat[[4]])) %>% 
              mutate_all(funs(unlist)) %>% 
              mutate(method = "average",
                     cluster = 1:10,
                     nr = "Solution",
                     aggregate = c(NA, unlist(aggregated_indizes[[17]][4,-1]))) )  


a1 = bench_results_df %>% 
  #filter(method!="average", method!="kmeans", method!="pam") %>%
  filter(method!="average") %>% 
  ggplot(aes(x=cluster, y=aggregate, col=method, shape=method)) +
  geom_line() +
  geom_point() +
  geom_point(size=3) +
  scale_color_grey(end = 0.5) +
  ylab("Aggregated Indices") +
  scale_x_continuous(breaks = seq(2,10,1)) +
  geom_hline(yintercept = 0) +
  theme_bw()
a2 = bench_results_df %>% 
  #filter(method!="average", method!="kmeans", method!="pam") %>% 
  filter(method!="average") %>% 
  ggplot(aes(x=cluster, y=avewithin, col=method, shape=method)) +
  geom_line() +
  geom_point() +
  geom_point(size=3) +
  scale_color_grey(end = 0.5) +
  ylab("Average Within Distance") +
  scale_x_continuous(breaks = seq(2,10,1)) +
  geom_hline(yintercept = 0) +
  theme_bw()
a3 = bench_results_df %>% 
  #filter(method!="average", method!="kmeans", method!="pam") %>% 
  filter(method!="average") %>% 
  ggplot(aes(x=cluster, y=pearsongamma, col=method, shape=method)) +
  geom_line() +
  geom_point() +
  geom_point(size=3) +
  scale_color_grey(end = 0.5) +
  ylab(expression("Pearson" * Gamma)) +
  scale_x_continuous(breaks = seq(2,10,1)) +
  geom_hline(yintercept = 0) +
  theme_bw()
a4 = bench_results_df %>% 
  #filter(method!="average", method!="kmeans", method!="pam") %>%
  filter(method!="average") %>% 
  ggplot(aes(x=cluster, y=asw, col=method, shape=method)) +
  geom_line() +
  geom_point(size=3) +
  scale_color_grey(end = 0.5) +
  ylab("Average Silhoutte Width") +
  scale_x_continuous(breaks = seq(2,10,1)) +
  geom_hline(yintercept = 0) +
  theme_bw()

ggarrange(a1,a2,a3, common.legend = T)  %>% 
  annotate_figure(top="Z-score calibration based on the same K")

# Stability ####
# 3 Cluster Solution
cboot_PAM_3_boot <- clusterboot(cluster_data,B=100, bootmethod=
                             c("boot"),clustermethod=claraCBI,
                           k=3, seed=15555)

cboot_PAM_3_boot
round(cboot_PAM_3_boot$bootmean, 3)

cboot_PAM_3_subset50 <- clusterboot(cluster_data,B=100, bootmethod=
                             c("subset"),clustermethod=claraCBI,
                           k=3, seed=15555, subtuning=nrow(cluster_data)*0.5)

cboot_PAM_3_subset50
round(cboot_PAM_3_subset50$subsetmean, 3)

cboot_PAM_3_subset25 <- clusterboot(cluster_data,B=100, bootmethod=
                             c("subset"),clustermethod=claraCBI,
                           k=3, seed=15555, subtuning=nrow(cluster_data)*0.25)

cboot_PAM_3_subset25
round(cboot_PAM_3_subset25$subsetmean, 3)

cboot_PAM_3_subset125 <- clusterboot(cluster_data,B=100, bootmethod=
                             c("subset"),clustermethod=claraCBI,
                           k=3, seed=15555, subtuning=nrow(cluster_data)*0.125)

cboot_PAM_3_subset125
round(cboot_PAM_3_subset125$subsetmean, 3)


# 4 Cluster Solution
cboot_PAM_4_boot <- clusterboot(cluster_data,B=100,bootmethod=
                             c("boot"),clustermethod=claraCBI,
                           k=4, seed=15555)

cboot_PAM_4_boot
round(cboot_PAM_4_boot$bootmean, 3)

cboot_PAM_4_subset50 <- clusterboot(cluster_data,B=100, bootmethod=
                                      c("subset"),clustermethod=claraCBI,
                                    k=4, seed=15555, subtuning=nrow(cluster_data)*0.5)

cboot_PAM_4_subset50
round(cboot_PAM_4_subset50$subsetmean, 3)

cboot_PAM_4_subset25 <- clusterboot(cluster_data,B=100, bootmethod=
                                      c("subset"),clustermethod=claraCBI,
                                    k=4, seed=15555, subtuning=nrow(cluster_data)*0.25)

cboot_PAM_4_subset25
round(cboot_PAM_4_subset25$subsetmean, 3)

cboot_PAM_4_subset125 <- clusterboot(cluster_data,B=100, bootmethod=
                                       c("subset"),clustermethod=claraCBI,
                                     k=4, seed=15555, subtuning=nrow(cluster_data)*0.125)

cboot_PAM_4_subset125
round(cboot_PAM_4_subset125$subsetmean, 3)

# 5 Cluster Solution
cboot_PAM_5_boot <- clusterboot(cluster_data,B=100,bootmethod=
                                  c("boot"),clustermethod=claraCBI,
                                k=5, seed=15555)

cboot_PAM_5_boot
round(cboot_PAM_5_boot$bootmean, 3)

cboot_PAM_5_subset50 <- clusterboot(cluster_data,B=100, bootmethod=
                                      c("subset"),clustermethod=claraCBI,
                                    k=5, seed=15555, subtuning=nrow(cluster_data)*0.5)

cboot_PAM_5_subset50
round(cboot_PAM_4_subset50$subsetmean, 3)

cboot_PAM_5_subset25 <- clusterboot(cluster_data,B=100, bootmethod=
                                      c("subset"),clustermethod=claraCBI,
                                    k=5, seed=15555, subtuning=nrow(cluster_data)*0.25)

cboot_PAM_5_subset25
round(cboot_PAM_5_subset25$subsetmean, 3)

cboot_PAM_5_subset125 <- clusterboot(cluster_data,B=100, bootmethod=
                                       c("subset"),clustermethod=claraCBI,
                                     k=5, seed=15555, subtuning=nrow(cluster_data)*0.125)

cboot_PAM_5_subset125
round(cboot_PAM_5_subset125$subsetmean, 3)



# 9 Cluster Solution
cboot_PAM_9_boot <- clusterboot(cluster_data,B=100,bootmethod=
                             c("boot"),clustermethod=claraCBI,
                           k=9, seed=15555)

cboot_PAM_9_boot
round(cboot_PAM_9_boot$bootmean, 3)


cboot_PAM_9_subset75 <- clusterboot(cluster_data,B=100, bootmethod=
                                      c("subset"),clustermethod=claraCBI,
                                    k=9, seed=15555, subtuning=nrow(cluster_data)*0.75)

cboot_PAM_9_subset75

round(cboot_PAM_9_subset50$subsetmean, 3)
cboot_PAM_9_subset50 <- clusterboot(cluster_data,B=100, bootmethod=
                                      c("subset"),clustermethod=claraCBI,
                                    k=9, seed=15555, subtuning=nrow(cluster_data)*0.5)

cboot_PAM_9_subset50
round(cboot_PAM_9_subset50$subsetmean, 3)

cboot_PAM_9_subset25 <- clusterboot(cluster_data,B=100, bootmethod=
                                      c("subset"),clustermethod=claraCBI,
                                    k=9, seed=15555, subtuning=nrow(cluster_data)*0.25)

cboot_PAM_9_subset25
round(cboot_PAM_9_subset25$subsetmean, 3)

cboot_PAM_9_subset125 <- clusterboot(cluster_data,B=100, bootmethod=
                                       c("subset"),clustermethod=claraCBI,
                                     k=9, seed=15555, subtuning=nrow(cluster_data)*0.125)

cboot_PAM_9_subset125
round(cboot_PAM_9_subset125$subsetmean, 3)



### Plot Results ####

FKM_2 = FKM(X = cluster_data, 2, maxit=1400, stand=0, RS=30, seed=1234)
FKM_3 = FKM(X = cluster_data, 3, maxit=1400, stand=0, RS=30, seed=1234)
FKM_4 = FKM(X = cluster_data, 4, maxit=1400, stand=0, RS=30, seed=1234)
FKM_9 = FKM(X = cluster_data, 9, maxit=2000, stand=0, RS=30, seed=1234)

PAM_2 = pamk(cluster_data, 2, seed=1234)
PAM_3 = pamk(cluster_data, 3, seed=1234)
PAM_4 = pamk(cluster_data, 4, seed=1234)
PAM_9 = pamk(cluster_data, 9, seed=1234)

kmeans_2 = kmeans(cluster_data, 2, nstart =20)
kmeans_3 = kmeans(cluster_data, 3, nstart =20)
kmeans_4 = kmeans(cluster_data, 4, nstart =20)
kmeans_9 = kmeans(cluster_data, 9, nstart =20, iter.max = 100)


# Principal Component Plot ####
principal_plot = function(cluster, label) {
  data.frame(predict(prcomp(cluster_data))[,1:2]) %>%  
    bind_cols(data.frame(Cluster = as.factor(cluster))) %>% 
    ggplot(aes(x=PC1, y=PC2, col=Cluster)) +
    geom_point(size=0.75) +
    ggtitle(label) +
    theme_bw()  +
    theme(legend.position = c("none")) 
}

principal_plot_FKM = function(cluster, label) {
  data.frame(predict(prcomp(cluster_data))[,1:2]) %>%  
    bind_cols(data.frame(Cluster = as.factor(cluster[,1]),
                         membership = cluster[,2])) %>% 
    ggplot(aes(x=PC1, y=PC2, col = Cluster, 
               alpha=membership)) +
    geom_point(size=0.75) +
    ggtitle(label) +
    theme_bw()  +
    theme(legend.position = c("none")) 
}


p1 = principal_plot_FKM(FKM_3$clus, "FCM 3")
p2 = principal_plot_FKM(FKM_4$clus, "FCM 4")
p3 = principal_plot_FKM(FKM_9$clus, "FCM 9")
p4 = principal_plot(PAM_3$pamobject$clustering, "PAM 3")
p5 = principal_plot(PAM_4$pamobject$clustering, "PAM 4")
p6 = principal_plot(PAM_9$pamobject$clustering, "PAM 9")
p7 = principal_plot(kmeans_3$cluster, "KMEANS 3")
p8 = principal_plot(kmeans_4$cluster, "KMEANS 4")
p9 = principal_plot(kmeans_9$cluster, "KMEANS 9")
ggarrange(p1,p2,p4,p5,p7,p8, ncol=3,nrow=3, common.legend = F)

ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9, ncol=3,nrow=3, common.legend = F)
# Explained Variability
round(summary(prcomp(cluster_data))$importance[3,2] * 100, 0)


# the cluster in the middle resembles the balanced cluster
p8 +
  theme(legend.position = c("bottom")) 


# Dimensions Plot ####
dim_plot = function(cluster, label) {
  dim_data_plot = dmx_data  %>%
    filter(classification_context == "Working Democracy" | classification_context == "Deficient Democracy") %>% 
    select(freedom = freedom_dim_index_context,
           equality = equality_dim_index_context,
           control = control_dim_index_context) %>% 
    bind_cols(data.frame(Cluster = as.factor(cluster))) %>% 
    mutate(Cluster = fct_reorder(Cluster, freedom)) 
  
  if (length(unique(dim_data_plot$Cluster)) == 2) {
    # dimlabels =  c("FeC", "FEC")
    } 
  if (length(unique(dim_data_plot$Cluster)) == 3) {
    levels(dim_data_plot$Cluster) =   c("Illiberal", "Inegalitarian", "Unaccountable")
    dim_data_plot$Cluster = fct_relevel(dim_data_plot$Cluster,  "Illiberal", "Inegalitarian", "Unaccountable")
  } 
  if (length(unique(dim_data_plot$Cluster)) == 4) {
    levels(dim_data_plot$Cluster) =  c("Illiberal", "Unaccountable", "Inegalitarian", "Balanced")
    dim_data_plot$Cluster = fct_relevel(dim_data_plot$Cluster,  "Illiberal", "Inegalitarian", "Unaccountable", "Balanced")
    
  }
  
  dim_data_plot %>% 
    pivot_longer(cols=c("freedom", "equality", "control")) %>% 
    mutate(name = fct_relevel(name, "freedom","equality","control")) %>% 
    ggplot(aes(x=Cluster, y=value, fill=name))+
    geom_boxplot() +
    stat_summary(fun.data = stat_box_count, geom = "text", hjust = 0.5, vjust = 0.6) + 
    # scale_x_discrete(labels = dimlabels) +
    scale_fill_grey(name="Dimensions", start = 0.4, end = 0.85) +
    theme_bw() +
    xlab("") +
    ylab("") +
    ggtitle(label) 
}
  
p11_dim = dim_plot(FKM_2$clus[,1], "FCM 2")
p1_dim = dim_plot(FKM_3$clus[,1], "FCM 3")
p2_dim = dim_plot(FKM_4$clus[,1], "FCM 4")
p3_dim = dim_plot(PAM_2$pamobject$clustering, "PAM 2")
p4_dim = dim_plot(PAM_3$pamobject$clustering, "PAM 3")
p5_dim = dim_plot(PAM_4$pamobject$clustering, "PAM 4")
p6_dim = dim_plot(kmeans_2$cluster, "KMEANS 2")
p7_dim = dim_plot(kmeans_3$cluster, "KMEANS 3")
p8_dim = dim_plot(kmeans_4$cluster, "KMEANS 4")
#p9_dim = dim_plot(kmeans_9$cluster, "KMEANS 9")
ggarrange(p1_dim,p2_dim,p4_dim,p5_dim,p7_dim,p8_dim,
          ncol=2,nrow=3, common.legend = TRUE, legend="bottom")


# Four Dim Plot
dim4_plot = function(cluster, label) {
  dim_data_plot = cluster_data  %>%
    select(freedom = freedom_uv,
           equality = equality_uv,
           control = control_uv,
           level = mean_dim_uv) %>%
    bind_cols(data.frame(Cluster = as.factor(cluster))) %>% 
    mutate(Cluster = fct_reorder(Cluster, level)) 
  
  if (length(unique(dim_data_plot$Cluster)) == 2) {
    # dimlabels =  c("FeC", "FEC")
    dimlabels =  c("Inegalitarian", "Balanced")
  } 
  if (length(unique(dim_data_plot$Cluster)) == 3) {
    # dimlabels =  c("fEC", "FeC", "FEC")    
    levels(dim_data_plot$Cluster) =   c("Illiberal", "Inegalitarian", "Unaccountable")
    dim_data_plot$Cluster = fct_relevel(dim_data_plot$Cluster,  "Illiberal", "Inegalitarian", "Unaccountable")
  } 
  if (length(unique(dim_data_plot$Cluster)) == 4) {
    # dimlabels =  c("fEC", "FEc", "FeC", "FEC")
    # dimlabels =  c("Illiberal", "Inegalitarian", "Unaccountable", "Balanced")
    levels(dim_data_plot$Cluster) =  c("Illiberal", "Inegalitarian", "Unaccountable", "Balanced")
    dim_data_plot$Cluster = fct_relevel(dim_data_plot$Cluster,  "Illiberal", "Inegalitarian", "Unaccountable", "Balanced")
  }
  
  
  
  dim_data_plot %>% 
    pivot_longer(cols=c("freedom", "equality", "control", "level")) %>% 
    mutate(name = fct_relevel(name, "freedom","equality","control")) %>% 
    ggplot(aes(x=Cluster, y=value, fill=name))+
    geom_boxplot() +
    stat_summary(fun.data = stat_box_count, geom = "text", hjust = 0.5, vjust = 2.1) + 
    # scale_x_discrete(labels = dimlabels) +
    scale_fill_grey(name="Dimensions", start = 0.4, end = 1) +
    theme_bw() +
    xlab("") +
    ylab("") +
    ylim(-0.1, 1) +
    ggtitle(label) 
}


p11_dim4 = dim4_plot(FKM_2$clus[,1], "FCM 2")
p1_dim4 = dim4_plot(FKM_3$clus[,1], "FCM 3")
p2_dim4 = dim4_plot(FKM_4$clus[,1], "FCM 4")
p3_dim4 = dim4_plot(PAM_2$pamobject$clustering, "PAM 2")
p4_dim4 = dim4_plot(PAM_3$pamobject$clustering, "PAM 3")
p5_dim4 = dim4_plot(PAM_4$pamobject$clustering, "PAM 4")
p6_dim4 = dim4_plot(kmeans_2$cluster, "KMEANS 2")
p7_dim4 = dim4_plot(kmeans_3$cluster, "KMEANS 3")
p8_dim4 = dim4_plot(kmeans_4$cluster, "KMEANS 4")
#p9_dim = dim_plot(kmeans_9$cluster, "KMEANS 9")
ggarrange(p1_dim4,p2_dim4,p4_dim4,p5_dim4,p7_dim4,p8_dim4,
          ncol=2,nrow=3, common.legend = TRUE, legend="bottom")


# Create and Save Cluster Dataset ####
complete_data_dimensions = dmx_data  %>%
  as.data.frame() %>% 
  filter(classification_context == "Working Democracy" | classification_context == "Deficient Democracy") %>% 
  rename(freedom = freedom_dim_index_context,
         equality = equality_dim_index_context,
         control = control_dim_index_context) %>%  
  bind_cols(data.frame(Cluster = as.factor(FKM_4$clus[,1]))) %>%  
  bind_cols(data.frame(mp = FKM_4$U)) %>% 
  mutate(Cluster = fct_recode(Cluster,
                              "Illiberal Democracy" = "2",
                              "Unaccountable Democracy" = "1",
                              "Inegalitarian Democracy" = "4",
                              "Balanced Democracy" = "3")) %>%
  rename(mp_Illiberal = mp.Clus.2,
         mp_Inegalitarian = mp.Clus.4,
         mp_Balanced = mp.Clus.3,
         mp_Unaccountable = mp.Clus.1) %>% 
  select_at(vars(country, year, regions, classification_context, freedom, equality, control, Cluster, mp_Illiberal, mp_Inegalitarian, mp_Balanced, mp_Unaccountable)) 

# write.csv(complete_data_dimensions, "Datasets/dmx_cluster_data_v3.csv")
complete_data_dimensions = read.csv("Datasets/dmx_cluster_data_v3.csv")


# Development Plot ####
plot_types_N = complete_data_dimensions %>%
  select(Cluster, year) %>%
  group_by(year) %>%
  summarise(n_total=n())

complete_data_cluster = data.frame(Cluster = rep(unique(complete_data_dimensions$Cluster), each=length(unique(complete_data_dimensions$year))),
                                   year = rep(unique(complete_data_dimensions$year), length(unique(complete_data_dimensions$Cluster))))

plot_types_yearly = complete_data_dimensions %>%
  select(Cluster, year) %>%
  group_by(year, Cluster) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  full_join(complete_data_cluster, by=c("year", "Cluster")) %>%
  arrange(year, Cluster) %>%
  mutate(
    n = ifelse(is.na(n) == T, 0, n)
  ) %>%
  left_join(plot_types_N, by="year") %>%
  mutate(percent = n/n_total)



dev1 = ggplot(plot_types_yearly, aes(x=year, y=n, fill=Cluster, shape=Cluster)) + 
  geom_area(stat="identity", col="black") + 
  theme_classic() +
  scale_x_continuous(breaks=seq(1900, 2020, 20)) + 
  theme(legend.position = "bottom", axis.text.x = element_text(angle=90, size=10), plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5)) + 
  xlab("") +
  coord_cartesian(expand=0) + 
  scale_fill_brewer(name="", type="qual", palette="Paired") + 
  #ggtitle("Temporal Distribution of Subtypes of Democracy", subtitle = "Count") + 
  ylab("Count") 

dev2 = ggplot(plot_types_yearly, aes(x=year, y=percent, fill=Cluster)) + 
  geom_area(stat="identity", col="black", size=0.8) + 
  theme_classic() +
  scale_x_continuous(breaks=seq(1900, 2000, 20), limits=c(1900, 2020)) + 
  theme(legend.position = "bottom", axis.text.x = element_text(angle=90, size=10), plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5)) + xlab("") +
  coord_cartesian(expand=0) + 
  scale_fill_brewer(name="", type="qual", palette="Paired") + 
  #ggtitle("Temporal Distribution of Subtypes of Democracy", subtitle = "Percent") + 
  scale_y_continuous(labels=percent, name="")


ggarrange(dev1, dev2, common.legend = T, nrow=1, ncol=2, legend="bottom") %>% 
  annotate_figure("Temporal Distribution of Subtypes of Democracy")

# Cross-Tabulation ####
library(pixiedust)
library(broom)

with(complete_data_dimensions, table(classification_context, Cluster)) %>% 
  tidy() %>% 
  group_by(classification_context) %>% 
  mutate(percent = n/sum(n),
         percent = round(percent,3) * 100,
         n = paste(n, " (", percent, ")", sep="")) %>% 
  select(-percent) %>% 
  pivot_wider(names_from = classification_context	, values_from = n) %>% 
  #arrange(classification_context) %>% 
  dust() %>% 
  sprinkle_colnames("Cluster", "Deficient \nDemocracy", "Working \nDemocracy") %>% 
  sprinkle_print_method("html")

# Regional Distribution ####

complete_data_dimensions %>% 
  group_by(regions, Cluster) %>% 
  summarise(no = n()) %>% 
  group_by(regions) %>% 
  mutate(perc = no/sum(no)) %>% 
  ungroup() %>% 
  mutate(regions = gsub(" ", "\n", regions)) %>%
  ggplot(aes(x=regions, y=perc, fill=Cluster)) +
  geom_bar(stat="identity") +
  scale_y_continuous(labels=scales::percent, name=NULL) +
  xlab("")  +
  theme_bw() +
  theme(axis.text.x = element_text(size=12))
 

# World Map ####

create_world_map= function(selected_year, label) {
  library(RColorBrewer)
  library(rworldmap)
  
  dmx_year = complete_data_dimensions %>% 
    filter(year==selected_year) %>% 
    dplyr::select(country, variable = Cluster) %>%
    mutate(country = as.character(country)) 
  
  prev_year_data = anti_join(complete_data_dimensions %>% 
                               filter(year==selected_year-1) %>% 
                               dplyr::select(country, variable = Cluster) %>%
                               mutate(country = as.character(country)), dmx_year, by="country")
  print(prev_year_data)
  
  dmx_year = dmx_year %>% 
    bind_rows(prev_year_data)
  
 
  
  dmx_year$country[dmx_year$country=="Burma/Myanmar"] = "Burma"
  dmx_year$country[dmx_year$country=="Republic of Vietnam"] = "Vietnam"
  dmx_year$country[dmx_year$country=="São Tomé and Príncipe"] = "Sao Tome and Principe"
  
  merged_map_data <- joinCountryData2Map(dmx_year,
                                         joinCode = "NAME",
                                         nameJoinColumn = "country",
                                         verbose = TRUE)
  cnt = as.character(merged_map_data$NAME[merged_map_data$NAME != "Antarctica"])
  cnt = as.character(cnt[cnt != "Greenland"])
  merged_map_data <- subset(merged_map_data, NAME  %in%  cnt)
  

  colourPalette <- brewer.pal(length(levels(dmx_year$variable)), 'RdYlBu')
  
  par(mai=c(0.6,0.1,0.6,0.1),xaxs="i",yaxs="i")
  mapParams = mapCountryData(merged_map_data, 
                             nameColumnToPlot="variable", 
                             colourPalette=colourPalette,
                             catMethod="categorical", 
                             addLegend = T,
                             borderCol= "black",
                             lwd=1, 
                             mapTitle = paste(label, selected_year),
                             missingCountryCol="lightgrey",
                             #mapRegion = "Europe"
  )
  do.call( addMapLegendBoxes, c(mapParams,title="Subtype"))
}
create_world_map(2019, "Regional Distribution of Subtypes of Democracy \n")
create_world_map(2018, "Regional Distribution of Subtypes of Democracy \n")
create_world_map(2017, "Regional Distribution of Subtypes of Democracy \n")
create_world_map(1960, "Regional Distribution of Subtypes of Democracy \n")




# 3D Tenery Plot ####

# Compute tetrahedron coordinates according to https://mathoverflow.net/a/184585
library(plot3D)
library(geometry)

simplex <- function(n) {
  qr.Q(qr(matrix(1, nrow=n)) ,complete = TRUE)[,-1]
}
tetra <- simplex(4)

df = complete_data_dimensions %>% 
  select_at(vars(starts_with("mp_"))) %>% 
  mutate_at(vars(starts_with("mp_")), funs(. + rnorm(length(.), 0, 0.02))) %>% 
  mutate(mysum = rowSums(.)) %>% 
  mutate_at(vars(starts_with("mp_")), funs(./mysum)) %>%
  select(-mysum) %>% 
  rename_at(vars(starts_with("mp_")), funs(gsub("mp_","",.))) %>% 
  #sample_n(500) %>% 
  as.matrix() 

# Convert barycentric coordinates (4D) to cartesian coordinates (3D)
df3D <- bary2cart(tetra, df)
labels = c("\n\n Illiberal\nDemocracy", "  Inegalitarian\n  Democracy", "  Balanced\n  Democracy", " Unaccountable\n    Democracy")


# Plot data

scatter3D(df3D[,1], df3D[,2], df3D[,3],
          xlim = range(tetra[,1]), ylim = range(tetra[,2]), zlim = range(tetra[,3]),
          col = "black", alpha = 0.075, pch = 16, box = FALSE, theta = 55, phi=5)
lines3D(tetra[c(1,2,3,4,1,3,1,2,4),1],
        tetra[c(1,2,3,4,1,3,1,2,4),2],
        tetra[c(1,2,3,4,1,3,1,2,4),3],
        col = "black", add = TRUE, lwd = 3, lty=2, alpha=0.4)
text3D(tetra[,1], tetra[,2], tetra[,3],
       labels, cex=1.2, add = TRUE)

