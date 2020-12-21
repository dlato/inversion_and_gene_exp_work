#test to see if coef of var are different between groups
#https://cran.r-project.org/web/packages/cvequality/vignettes/how_to_test_CVs.html

library(cvequality)
library(dplyr)
library(plyr)

#read in dat
g_dat <- read.table("Stats/inversions_gene_exp_wtest_data.csv", sep = "\t")
summary(g_dat)

#CV test between inver and non-inver (all blocks)
# null = no difference in CV between groups
#Feltz and Miller's test (1996)
with(g_dat, 
     asymptotic_test(norm_exp,
                     inversion))
#Krishnamoorthy and Lee's modified signed-likelihood ratio test (2014)
with(g_dat, 
     mslr_test(nr = 1e4, 
               norm_exp,
               inversion))

#CV test between rev_comp and non-rev_comp (all blocks)
# null = no difference in CV between groups
#Feltz and Miller's test (1996)
with(g_dat, 
     asymptotic_test(norm_exp,
                     rev_comp))
#Krishnamoorthy and Lee's modified signed-likelihood ratio test (2014)
with(g_dat, 
     mslr_test(nr = 1e4, 
               norm_exp,
               rev_comp))

#CV test between SIG inver and non-SIG inver
sig_d <- g_dat %>% filter(!is.na(block_w_pvalue)) %>%
                  mutate(sig = ifelse(block_w_pvalue <= 0.05, "yes", "no"))
# null = no difference in CV between groups
#Feltz and Miller's test (1996)
with(sig_d, 
     asymptotic_test(norm_exp,
                     sig))
#Krishnamoorthy and Lee's modified signed-likelihood ratio test (2014)
with(sig_d, 
     mslr_test(nr = 1e4, 
               norm_exp,
               sig))

#CV test between rev_comp and non-rev_comp (ATCC blocks)
atcc_d <- g_dat %>% filter(strain == "ATCC")
# null = no difference in CV between groups
#Feltz and Miller's test (1996)
with(atcc_d, 
     asymptotic_test(norm_exp,
                     rev_comp))
#Krishnamoorthy and Lee's modified signed-likelihood ratio test (2014)
with(atcc_d, 
     mslr_test(nr = 1e4, 
               norm_exp,
               rev_comp))


# exp dist ALL BLOCKS
#format data
g_dat$inversion <- as.factor(g_dat$inversion)
g_mu <- ddply(g_dat, "inversion", summarise, grp.mean=mean(norm_exp))
p<- (ggplot(g_dat, aes(x=norm_exp, fill=inversion))
     + geom_density(alpha=0.4)
# Add mean lines
     + geom_vline(data=g_mu, aes(xintercept=grp.mean, color=inversion),
             linetype="dashed")
#     + scale_x_continuous(trans = 'log10')
   + labs(title= "Distribution of Gene Expression",x = "Gene Expression (CPM)", y = "Density") 
   + scale_x_continuous(trans = 'log10', labels = function(x) ifelse(x ==0, "0", x))
   + scale_y_continuous(labels = function(x) ifelse(x ==0,"0", x))
   + scale_color_manual(values = c("#788AA3","#2E294E"), labels = c("Inversion","Non-Inversion"))
   + scale_fill_manual(values = c("#788AA3","#2E294E"), labels = c("Inversion","Non-Inversion"))

)
pdf("Stats/all_exp_dist.pdf")
p
dev.off()

# exp dist ATCC
#format data
atcc_d$rev_comp <- as.factor(atcc_d$rev_comp)
summary(atcc_d)
a_mu <- ddply(atcc_d, "rev_comp", summarise, grp.mean=mean(norm_exp))
ecolT <- substitute(tit~italic(ecoli)~strain, list(tit = "Distribution of Gene Expression in", ecoli="E.coli",strain="ATCC 25922"))
p<- (ggplot(atcc_d, aes(x=norm_exp, fill=rev_comp))
     + geom_density(alpha=0.4)
# Add mean lines
     + geom_vline(data=a_mu, aes(xintercept=grp.mean, color=rev_comp),
             linetype="dashed")
#     + scale_x_continuous(trans = 'log10')
   + labs(title= ecolT,x = "Gene Expression (CPM)", y = "Density") 
   + scale_x_continuous(trans = 'log10', labels = function(x) ifelse(x ==0, "0", x))
   + scale_y_continuous(labels = function(x) ifelse(x ==0,"0", x))
   + scale_color_manual(values = c("#788AA3","#2E294E"), labels = c("Inversion","Non-Inversion"))
   + scale_fill_manual(values = c("#788AA3","#2E294E"), labels = c("Inversion","Non-Inversion"))

)
pdf("Stats/ATCC_exp_dist.pdf")
p
dev.off()
