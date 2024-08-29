library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(latex2exp)
library(egg)

pcs_oi <- c(0, 3, 5, 10, 20, 40)

gls_thin <- function(z_df, md) {
  i <- 2
  while (i < (nrow(z_df))) {
    if (!((abs(z_df$r2[i] - z_df$r2[i+1]) > md) & (abs(z_df$r2[i] - z_df$r2[i-1]) > md))) {
      z_df <- z_df[-i, ]
      i <- i - 1
    } 
    i <- i + 1
  }
  return(z_df)
  
}


get_breaks <- function(x, n_breaks = 4) {

  breaks <- seq(from = min(x), to = max(x), length.out=n_breaks)        
  
  thresh_breaks <- c()
  
  while(length(breaks) > 0) {
    break_val <- which.min(abs(x-breaks[1]))
    thresh_breaks <- c(thresh_breaks, x[break_val])
    x <- x[-break_val]
    breaks <- breaks[-1]
  }
  
  return(thresh_breaks)
  
}

get_breaks <- function(x, n_breaks = 4) {

  breaks <- seq(from = min(x), to = max(x), length.out=n_breaks)        
  
  thresh_breaks <- c()
  
  while(length(breaks) > 0) {
    break_val <- which.min(abs(x-breaks[1]))
    thresh_breaks <- c(thresh_breaks, x[break_val])
    x <- x[-break_val]
    breaks <- breaks[-1]
  }
  
  return(thresh_breaks)
  
}

format_thresh_labels <- function(thresh) {
  orig <- 10^(-1/thresh)
  
  if (orig >= 0.0001) {
    return(orig)
  } else {
    paste0("1e-", 1/thresh)
  }
}

main_hairpin <- function(full_df, as, pop, pcs = unique(full_df$pc_num)) {
    # produce a standard hairpin plot of all the phenotypes in full_df

    full_df <- full_df %>% filter(pc_num %in% pcs)
  
    if (as == "as") {
        title <- paste0("Hairpin, in ", pop, " correcting for assessment center")
    } else {
        title <- paste0("Hairpin, in ", pop, " without correcting for assessment center")
    }

    p <- ggplot(data=full_df, aes(y = theta_eo, x = r2, col=phenotype, size=1/-log10(threshold))) + 
        geom_point() + 
        geom_line(aes(linetype=factor(pc_num), col = phenotype), linewidth=0.25) + 
        scale_size_binned(breaks = c(1, 1/2, 1/5, 1/10, 1/30, 1/60), range = c(0.4, 3.5), labels = function(x) {return(paste0("1e-", 1/x))} ) +
        ggtitle(title) + 
        labs(linetype = "PC", size = "Threshold") + 
        ylab(TeX("$\\theta_{eo}$")) +
        xlab(TeX("$\\r^2$")) +
        ylim(min(c(full_df$theta_eo, 0)), NA) +
        xlim(min(c(full_df$r2, 0)), NA) +
        theme_classic() +
        scale_colour_brewer(palette = "Dark2") 

    return(p)
}

main_hairpin_flip <- function(full_df, as, pop, pcs = unique(z_df$pc_num)) {
    # flipped axis of main_hairpin
  
    full_df <- full_df %>% filter(pc_num %in% pcs)
  
    if (as == "as") {
        title <- paste0("Hairpin, in ", pop, " correcting for assessment center")
    } else {
        title <- paste0("Hairpin, in ", pop, " without correcting for assessment center")
    }

    p <- ggplot(data=full_df, aes(y = r2, x = theta_eo, col=phenotype, size=1/-log10(threshold))) + 
        geom_point() + 
        geom_line(aes(linetype=factor(pc_num), col = phenotype), linewidth=0.25) + 
        scale_size_binned(breaks = c(1, 1/2, 1/5, 1/10, 1/30, 1/60), range = c(0.4, 3.5), labels = function(x) {return(paste0("1e-", 1/x))} ) +
        ggtitle(title) + 
        labs(linetype = "PC", size = "Threshold") + 
        ylab(TeX("$\\r^2$")) +
        xlab(TeX("$\\theta_{eo}$")) +
        ylim(min(c(full_df$r2, 0)), NA) +
        xlim(min(c(full_df$theta_eo, 0)), NA) +
        theme_classic() +
        scale_colour_brewer(palette = "Dark2") 

  return(p)
}
  
indiv_hairpin <- function(p_df, as, phen, pop, pcs = unique(p_df$pc_num), n_breaks = 5) {
    # produces a hairpin plot for an individual phenotype

    p_df <- p_df %>% filter(pc_num %in% pcs)

    if (as == "as") {
        title <- paste0("Hairpin for ", phen, ", in ", pop, " correcting for assessment center")
    } else {
        title <- paste0("Hairpin for ", phen, ", in ", pop, " without correcting for assessment center")
    }

    bs <- unique(p_df$threshold)  

    if (1 %in% bs) {
        bs <- bs[-which(bs == 1)]
        bs <- 1/-log10(unique(bs))

        brks <- get_breaks(bs, n_breaks)

        labs <- as.numeric(sapply(brks, FUN = (format_thresh_labels)))
        labs[which.max(as.numeric(labs))] = 1
        labs <- labs[order(labs)]
    } else {
        bs <- 1/-log10(unique(p_df$threshold))

        brks <- get_breaks(bs, n_breaks)

        labs <- as.numeric(sapply(brks, FUN = (format_thresh_labels)))
        labs <- labs[order(labs)]
    }


    p <- ggplot(data=p_df, aes(y = theta_eo, x = r2, size=case_when(threshold == 1 ~ 1, .default = 1/-log10(threshold)))) + 
        geom_point(colour="#E66100") + 
        geom_line(aes(linetype=factor(pc_num)), linewidth=0.25, colour = "#5D3A9B") + 
        scale_size_binned(breaks = brks, range = c(1, 3), labels = labs) +
        ggtitle(title) + 
        labs(linetype = "PC", size = "Threshold") + 
        ylab(TeX("$\\theta_{eo}$")) +
        xlab(TeX("$\\r^2$")) +
        ylim(min(c(p_df$theta_eo, 0)), NA) +
        xlim(min(c(p_df$r2, 0)), NA) +
        theme_classic() +
        scale_colour_brewer(palette = "Dark2")

    return(p)
} 

indiv_hairpin_flip <- function(p_df, as, phen, pop, pcs = unique(p_df$pc_num), n_breaks = 5) {
    # flipped axis for indiv_hairpin

    p_df <- p_df %>% filter(pc_num %in% pcs)

    if (as == "as") {
        title <- paste0("Hairpin for ", phen, ", in ", pop, " correcting for assessment center")
    } else {
        title <- paste0("Hairpin for ", phen, ", in ", pop, " without correcting for assessment center")
    }

    bs <- unique(p_df$threshold)  

    if (1 %in% bs) {
        bs <- bs[-which(bs == 1)]
        bs <- 1/-log10(unique(bs))

        brks <- get_breaks(bs, n_breaks)

        labs <- as.numeric(sapply(brks, FUN = (format_thresh_labels)))
        labs[which.max(as.numeric(labs))] = 1
        labs <- labs[order(labs)]
    } else {
        bs <- 1/-log10(unique(p_df$threshold))

        brks <- get_breaks(bs, n_breaks)

        labs <- as.numeric(sapply(brks, FUN = (format_thresh_labels)))
        labs <- labs[order(labs)]
    }

    p <- ggplot(data=p_df, aes(y = r2, x = theta_eo, size=case_when(threshold == 1 ~ 1, .default = 1/-log10(threshold)))) + 
        geom_point(colour="#E66100") + 
        geom_line(aes(linetype=factor(pc_num)), linewidth=0.25, colour = "#5D3A9B") + 
        scale_size_binned(breaks = brks, range = c(1, 3), labels = labs) +
        ggtitle(title) + 
        labs(linetype = "PC", size = "Threshold") + 
        ylab(TeX("$\\r^2$")) +
        xlab(TeX("$\\theta_{eo}$")) +
        ylim(min(c(p_df$r2, 0)), NA) +        
        xlim(min(c(p_df$theta_eo, 0)), NA) +
        theme_classic() +
        scale_colour_brewer(palette = "Dark2")

    return(p)
} 

cols <- brewer.pal(n=10, name='Paired')
pcs <- c("0", "1", "2", "3", "4", "5", "10", "20", "30", "40")

map_pc <- function(pc) {
  return(cols[which(pcs==pc)])
}

indiv_gls_p <- function(z_df, phen, pc_n, as = "as", pop = "WBRT", md = 0) {
    # produces a plot of an individual GLS line with points labeled by threshold

    if (as == "as") {
        title <- paste0("Hairpin for ", phen, ", PC ", pc_n, " (GLS), in ", pop, " correcting for assessment center. Labeled by thresh")
    } else {
        title <- paste0("Hairpin for ", phen, ", PC ", pc_n, " (GLS), in ", pop, "w/out correcting for assessment center. Labeled by thresh")
    }

    z_df <- z_df %>% filter(pc_num == pc_n)

    z_df <- gls_thin(z_df, md)

    p <- ggplot(z_df, aes(x=r2, y=theta_eo, colour = on_line, label = threshold)) + 
        geom_point() +
        geom_line(colour = "black") +
        geom_abline(aes(intercept = 0, slope = beta), color = "blue", linewidth = .5) +
        geom_text(hjust = 1, vjust = 1.2, show.legend = FALSE) +
        ggtitle(title) +
        xlim(min(z_df$r2, 0), NA) + 
        ylim(min(z_df$theta_eo, 0), NA) +
        ylab(TeX("$\\theta_{eo}$")) +
        xlab(TeX("$\\r^2$")) +
        theme_classic() +
        scale_color_manual(values = c("on" = "black", "off" = "red"))

    return(p)
}

indiv_gls_sd <- function(z_df, phen, pc_n, as="as", pop = "WBRT", md = 0) {
    
    # produces a plot of an individual GLS line with points labeled by sd
    
    if (as == "as") {
        title <- paste0("Hairpin for ", phen, ", PC ", pc_n, " (GLS), in ", pop, " correcting for assessment center. Labeled by sd")
    } else {
        title <- paste0("Hairpin for ", phen, ", PC ", pc_n, " (GLS), in ", pop, "w/out correcting for assessment center. Labeled by sd")
    }
  
    z_df <- z_df %>% filter(pc_num == pc_n)
  
    z_df <- gls_thin(z_df, md)
  
    z_df$sd_eo <- format(z_df$sd_eo, scientific=TRUE, digits=3)

    p <- ggplot(z_df, aes(x=r2, y=theta_eo, colour = on_line, label = threshold)) + 
        geom_point() +
        geom_line(colour = "black") +
        geom_abline(aes(intercept = 0, slope = beta), color = "blue", linewidth = .5) +
        geom_text(hjust = 0, vjust = 4, show.legend = FALSE, size = 2.3) +
        ggtitle(title) +
        xlim(min(z_df$r2, 0), NA) + 
        ylim(min(z_df$theta_eo - as.numeric(max(z_df$sd_teo)), 0), max(z_df$theta_eo + as.numeric(max(z_df$sd_teo)))) +
        geom_errorbar(aes(ymin=as.numeric(theta_eo)-as.numeric(sd_teo), ymax=as.numeric(theta_eo)+as.numeric(sd_teo))) +
        ylab(TeX("$\\theta_{eo}$")) +
        xlab(TeX("$\\r^2$")) +
        theme_classic() +
        scale_color_manual(values = c("on" = "black", "off" = "red"))     

  return(p)
}

indiv_gls_CI <- function(z_df, phen, pc_n, as="as", pop = "WBRT", md = 0) {
    
    # produces a plot of an individual GLS line with points labeled by 95% CI

    if (as == "as") {
        title <- paste0("Hairpin for ", phen, ", PC ", pc_n, " (GLS), in ", pop, " correcting for assessment center. Labeled by 95% CI")
    } else {
        title <- paste0("Hairpin for ", phen, ", PC ", pc_n, " (GLS), in ", pop, "w/out correcting for assessment center. Labeled by 95% CI")
    }
  
    z_df <- z_df %>% filter(pc_num == pc_n)
  
    z_df <- gls_thin(z_df, md)
  
    z_df$sd_eo <- format(z_df$sd_eo, scientific=TRUE, digits=3)

    p <- ggplot(z_df, aes(x=r2, y=theta_eo, colour = on_line, label = threshold)) + 
        geom_point() +
        geom_line(colour = "black") +
        geom_abline(aes(intercept = 0, slope = beta), color = "blue", linewidth = .5) +
        geom_text(hjust = 0, vjust = 4, show.legend = FALSE, size = 2.3) +
        ggtitle(title) +
        xlim(min(z_df$r2, 0), NA) + 
        ylim(min(z_df$theta_eo - 1.96*as.numeric(max(z_df$sd_teo)), 0), max(z_df$theta_eo + 1.96*as.numeric(max(z_df$sd_teo)))) +
        ylab(TeX("$\\theta_{eo}$")) +
        xlab(TeX("$\\r^2$")) +
        geom_errorbar(aes(ymin=as.numeric(theta_eo)-1.96*as.numeric(sd_teo), ymax=as.numeric(theta_eo)+1.96*as.numeric(sd_teo))) +
        theme_classic() +
        scale_color_manual(values = c("on" = "black", "off" = "red"))     

  return(p)
}

all_gls_phen <- function(z_df, phen, as = "as", pop = "WBRT", pcs = unique(z_df$pc_num)) {
    
    # produces a plot of all the GLS lines for a given phenotype

    if (as == "as") {
        title <- paste0("Hairpin for ", phen, " (GLS), in ", pop, " correcting for assessment center")
    } else {
        title <- paste0("Hairpin for ", phen, " (GLS), in ", pop, "w/out correcting for assessment center")
    }
  
    z_df <- z_df %>% filter(pc_num %in% pcs)
   
    p <- ggplot(z_df, aes(x=r2, y=theta_eo, colour = factor(pc_num) )) + 
        geom_point(aes(pch = on_line, alpha = on_line, size = on_line)) +
        geom_line(aes(linetype = "solid", colour = factor(pc_num))) +
        geom_abline(aes(linetype = "longdash", colour=factor(pc_num), slope = beta, intercept = 0), linewidth = 0.25) + 
        labs(colour = "pc", size = "on_line", pch = "on_line") +
        guides(alpha="none") + 
        ggtitle(title) + 
        xlim(min(z_df$r2, 0), NA) + 
        ylim(min(z_df$theta_eo, 0), NA) + 
        ylab(TeX("$\\theta_{eo}$")) +
        xlab(TeX("$\\r^2$")) +
        theme_classic() +
        scale_colour_brewer(palette = "Paired") +
        scale_shape_manual(values = c("on" = 1, "off" = 16)) +
        scale_linetype_discrete(labels=c('GLS', 'Data')) +
        scale_color_hue(l=40, c=100) +
        scale_alpha_manual(values = c("on" = 1, "off" = 1)) +
        scale_size_manual(values = c("on" = 1.5, "off" = 3))
    
    return(p)
}


all_gls_phen_flip <- function(z_df, phen, as = "as", pop = "WBRT", pcs = unique(z_df$pc_num)) {
    
    # flipped axis for all_gls_phen

    if (as == "as") {
        title <- paste0("Hairpin for ", phen, "(GLS), in ", pop, " correcting for assessment center")
    } else {
        title <- paste0("Hairpin for ", phen, "(GLS), in ", pop, "w/out correcting for assessment center")
    }
  
    z_df <- z_df %>% filter(pc_num %in% pcs)
   
    p <- ggplot(z_df, aes(x=theta_eo, y=r2, colour = factor(pc_num) )) + 
        geom_point(aes(pch = on_line, alpha = on_line, size = on_line)) +
        geom_line(aes(linetype = "solid", colour = factor(pc_num))) +
        geom_abline(aes(linetype = "longdash", colour=factor(pc_num), slope = beta, intercept = 0), linewidth = 0.25) + 
        labs(colour = "pc", size = "on_line", pch = "on_line") +
        guides(alpha="none") + 
        ggtitle(title) + 
        xlim(min(z_df$theta_eo, 0), NA) + 
        ylim(min(z_df$r2, 0), NA) + 
        xlab(TeX("$\\theta_{eo}$")) +
        ylab(TeX("$\\r^2$")) +
        theme_classic() +
        scale_colour_brewer(palette = "Paired") +
        scale_shape_manual(values = c("on" = 1, "off" = 16)) +
        scale_linetype_discrete(labels=c('GLS', 'Data')) +
        scale_color_hue(l=40, c=100) +
        scale_alpha_manual(values = c("on" = 1, "off" = 1)) +
        scale_size_manual(values = c("on" = 1.5, "off" = 3))
    
    return(p)
}

all_gls_facet <- function(z_df, phen, as = "as", pop = "WBRT", pcs = unique(z_df$pc_num)) {
  
  if (as == "as") {
    title <- paste0("Hairpin for ", phen, " (GLS), in ", pop, " correcting for assessment center")
  } else {
    title <- paste0("Hairpin for ", phen, " (GLS), in ", pop, "w/out correcting for assessment center")
  }
  
  z_df <- z_df %>% filter(pc_num %in% pcs)

  p <- ggplot(z_df, aes(x=r2, y=theta_eo, colour = on_line)) + 
    geom_point() + 
    geom_line(colour = "black") +
    geom_abline(aes(slope = beta, intercept = 0), linewidth = 0.25, col = "blue") + 
    labs(pch = "pc", linetype = "pc") + 
    ggtitle(title) + 
    xlim(min(z_df$r2, 0), NA) + 
    ylim(min(z_df$theta_eo - max(z_df$sd_teo), 0), NA) +
    ylab(TeX("$\\theta_{eo}$")) +
    xlab(TeX("$\\r^2$")) +
    geom_errorbar(aes(ymin=as.numeric(theta_eo)-as.numeric(sd_teo), ymax=as.numeric(theta_eo)+as.numeric(sd_teo)), linewidth = 0.25) +
    theme_classic() +
    scale_color_manual(values = c("on" = "black", "off" = "red"))
    
  
  return(p + facet_wrap(vars(pc_num), ncol=3))
  
}


all_all_gls_phen <- function(zlist, as = "as", nrow = 2, ncol = 2) {
  p_df <- zlist[[1]]
  
  for (i in 2:length(zlist)) {
    p_df <- rbind(p_df, zlist[[i]])
  }
  
  p <- ggplot(p_df, aes(x=r2, y=theta_eo, colour = factor(pc_num) )) + 
    geom_point(aes(pch = on_line, alpha = on_line, size = on_line)) +
    geom_line(aes(linetype = "solid", colour = factor(pc_num))) +
    labs(colour = "pc", size = "on_line", pch = "on_line") + 
    geom_abline(aes(linetype = "longdash", colour=factor(pc_num), slope = beta, intercept = 0), linewidth = 0.25) + 
    xlim(min(p_df$r2, 0), NA) + 
    ylim(min(p_df$theta_eo, 0), NA) + 
    scale_shape_manual(values = c("on" = 1, "off" = 16)) +
    scale_linetype_discrete(labels=c('GLS', 'Data')) +
    scale_color_hue(l=40, c=100) +
    scale_alpha_manual(values = c("on" = 1, "off" = 1)) +
    scale_size_manual(values = c("on" = 1.5, "off" = 3)) +
    guides(alpha="none") + 
    theme_classic() +
    ylab(TeX("$\\theta_{eo}$")) +
    xlab(TeX("$\\r^2$")) + 
    scale_colour_brewer(palette = "Paired")
  
  return(p + facet_wrap(vars(phenotype), ncol=ncol, nrow = 2, scales = "free"))
}

all_all_gls_phen_inv <- function(zlist, as = "as", nrow = 2, ncol = 2) {
  p_df <- zlist[[1]]
  
  for (i in 2:length(zlist)) {
    p_df <- rbind(p_df, zlist[[i]])
  }
  
  p <- ggplot(p_df, aes(y=r2, x=theta_eo, colour = factor(pc_num) )) + 
    geom_point(aes(pch = on_line, alpha = on_line, size = on_line)) +
    geom_line(aes(linetype = "solid", colour = factor(pc_num))) +
    labs(colour = "pc", size = "on_line", pch = "on_line") + 
    geom_abline(aes(linetype = "longdash", colour=factor(pc_num), slope = 1/beta, intercept = 0), linewidth = 0.25) + 
    xlim(min(p_df$theta_eo, 0), NA) + 
    ylim(min(p_df$r2, 0), NA) + 
    scale_shape_manual(values = c("on" = 1, "off" = 16)) +
    scale_linetype_discrete(labels=c('GLS', 'Data')) +
    scale_color_hue(l=40, c=100) +
    scale_alpha_manual(values = c("on" = 1, "off" = 1)) +
    scale_size_manual(values = c("on" = 1.5, "off" = 3)) +
    guides(alpha="none") + 
    theme_classic() +
    xlab(TeX("$\\theta_{eo}$")) +
    ylab(TeX("$\\r^2$")) +
    scale_colour_brewer(palette = "Paired")
  
  return(p + facet_wrap(vars(phenotype), ncol=ncol, nrow = 2, scales = "free"))
}

mains	<- c( 'True Betas', 'Beta Estimates, Unadjusted', 'Beta Estimates, 1 PC Adjusted', 'Beta Estimates, 10 PC Adjusted' , 'Beta Estimates, True Pop Adjusted' )

sim_plot1 <- function(all_corves, gens, kk) {
  
  corves	<- apply( all_corves[,,,,,1], 1:4, mean, na.rm=T )
  
  corves_df_ve <- data.frame(corves['ve' ,kk,,])
  corves_df_ve$SNPS <- rownames(corves_df_ve)
  corves_df_ve <- melt(corves_df_ve)
  
  corves_df_ve <- corves_df_ve %>% rename(r2 = value)
  corves_df_ve <- corves_df_ve %>% rename(Fst = variable)
  
  corves_df_cor <- data.frame(corves['cor' ,kk,,])
  corves_df_cor$SNPS <- rownames(corves_df_cor)
  corves_df_cor <- melt(corves_df_cor)
  
  corves_df_cor <- corves_df_cor %>% rename(theta_eo = value)
  corves_df_cor <- corves_df_cor %>% rename(Fst = variable)
  
  full_df <- merge(corves_df_ve, corves_df_cor, by = c("SNPS", "Fst"))
  full_df <- full_df %>% mutate(gen = 1)
  
  rm(corves_df_ve)
  rm(corves_df_cor)
  
  for (i in 2:length(gens)) {
    corves	<- apply( all_corves[,,,,,i], 1:4, mean, na.rm=T )
    
    corves_df_ve <- data.frame(corves['ve' ,kk,,])
    corves_df_ve$SNPS <- rownames(corves_df_ve)
    corves_df_ve <- melt(corves_df_ve)
    
    corves_df_ve <- corves_df_ve %>% rename(r2 = value)
    corves_df_ve <- corves_df_ve %>% rename(Fst = variable)
    
    corves_df_cor <- data.frame(corves['cor' ,kk,,])
    corves_df_cor$SNPS <- rownames(corves_df_cor)
    corves_df_cor <- melt(corves_df_cor)
    
    corves_df_cor <- corves_df_cor %>% rename(theta_eo = value)
    corves_df_cor <- corves_df_cor %>% rename(Fst = variable)
    
    full_df_c <- merge(corves_df_ve, corves_df_cor, by = c("SNPS", "Fst"))
    full_df_c <- full_df_c %>% mutate(gen = i)
    
    rm(corves_df_ve)
    rm(corves_df_cor)
    
    full_df <- rbind(full_df, full_df_c)
    rm(full_df_c)
  }
  
  full_df$SNPS <- as.numeric(full_df$SNPS)
  
  if (kk == 1) {
    xlim1 <- min(0, full_df$r2)
    xlim2 <- 0.25
    ylim1 <- min(-0.02, full_df$theta_eo)
    ylim2 <- 0.05
  } else if (kk == 2) {
    xlim1 <- min(0, full_df$r2)
    xlim2 <- 0.04
    ylim1 <- min(-0.01, full_df$theta_eo)
    ylim2 <- 1
  } else if (kk == 3) {
    xlim1 <- min(0, full_df$r2)
    xlim2 <- 0.04
    ylim1 <- min(-0.01, full_df$theta_eo)
    ylim2 <- 1
  } else if (kk == 4) {
    xlim1 <- min(0, full_df$r2)
    xlim2 <- 0.04
    ylim1 <- min(-0.01, full_df$theta_eo)
    ylim2 <- 0.2
  } else {
    xlim1 <- min(0, full_df$r2)
    xlim2 <- 0.04
    ylim1 <- min(-0.01, full_df$theta_eo)
    ylim2 <- 0.2
  }
  
  p <- ggplot(data = full_df, aes(x=r2, y=theta_eo)) +
        geom_point(aes(size=factor(SNPS), colour = Fst)) +
        geom_line(aes(colour = Fst)) +
        ggtitle(paste0(mains[kk], " after ", length(gens), " generations of AM")) +
        scale_size_discrete(name = "# SNPs in PRS") +
        scale_colour_discrete(labels = c('1e-4', '0.001', '0.01', '0.1')) +
        theme_classic() +
        ylab(TeX("$\\theta_{eo}$")) +
        xlab(TeX("$\\r^2$")) +
        xlim(xlim1, xlim2) +
        ylim(ylim1, ylim2) + 
        facet_wrap(vars(gen), ncol=3) +
        scale_colour_brewer(palette = "Paired")
  
  return(p)
}


indiv_gls_sd_sim <- function(z_df, phen, pc_n, pc_lab, ngen) {
  
  title <- paste0("Hairpin for Fst=", phen, ",", pc_lab, "(GLS). Labeled by #SNPS")
  
  z_df <- z_df %>% filter(pc_num == pc_n)
  
  z_df$sd_eo <- format(z_df$sd_eo, scientific=TRUE, digits=3)
  beta <- unique(z_df$beta)[1]
  p <- ggplot(z_df, aes(x=r2, y=theta_eo, colour = on_line, label = threshold)) + 
    geom_line(colour = "black") +
    geom_point() +
    geom_text(hjust = 0, vjust = 4, show.legend = FALSE, size = 2.3) +
    theme_bw()  +
    labs(title = title, y = "theta E/O", x = "R2")  +
    scale_color_manual(values = c("on" = "black", "off" = "red")) + 
    geom_abline(intercept = 0, slope = beta, color = "blue", linewidth = .5) +
    geom_errorbar(aes(ymin=as.numeric(theta_eo)-as.numeric(sd_eo), ymax=as.numeric(theta_eo)+as.numeric(sd_eo))) +
    xlim(0, NA) + 
    ylim(0, NA) +
    theme_classic() +
    ylab(TeX("$\\theta_{eo}$")) +
    xlab(TeX("$\\r^2$"))
  
  return(p)
}


all_gls_phen_sim <- function(z_df, phen, ngen) {
  
  title <- paste0('Hairpin for Fst=', phen, ' (GLS). Labeled by #SNPS')
  
  p <- ggplot(z_df, aes(x=r2, y=theta_eo, colour = factor(pc_num) )) + 
    geom_point(aes(pch = on_line, alpha = on_line, size = on_line)) +
    geom_line(aes(linetype = "solid", colour = factor(pc_num))) +
    labs(colour = "pc", size = "on_line", pch = "on_line") + 
    geom_abline(aes(linetype = "longdash", colour=factor(pc_num), slope = 1/beta, intercept = 0), linewidth = 0.25) + 
    ggtitle(title) + 
    xlim(0, NA) + 
    ylim(0, NA) + 
    scale_shape_manual(values = c("on" = 1, "off" = 16)) +
    scale_linetype_discrete(labels=c('GLS', 'Data')) +
    scale_color_hue(l=40, c=100) +
    scale_alpha_manual(values = c("on" = 1, "off" = 1)) +
    scale_size_manual(values = c("on" = 1.5, "off" = 3)) +
    guides(alpha="none") + 
    theme_classic() +
    ylab(TeX("$\\theta_{eo}$")) +
    xlab(TeX("$\\r^2$")) +
    scale_colour_brewer(palette = "Paired")
  
  return(p)
}

indiv_dist_sd <- function(z_df, phen, pc_n, as="as") {
  
  if (as == "as") {
    title <- paste0("Hairpin for ", phen, ", PC ", pc_n, "(distance), w/ assessment center. Labeled by thresh")
  } else {
    title <- paste0("Hairpin for ", phen, ", PC ", pc_n, "(distance), w/out assessment center. Labeled by thresh")
  }
  
  z_df <- z_df %>% filter(pc_num == pc_n)
  
  p <- ggplot(z_df, aes(x=r2, y=mean_dist, colour = on_line, label = threshold)) + 
    geom_line(colour = "black") +
    geom_point() +
    geom_text(hjust = 0, vjust = 4, show.legend = FALSE, size = 2.3) +
    theme_bw()  +
    labs(title = title, y = "mean distance", x = "R2")  +
    scale_color_manual(values = c("on" = "black", "off" = "red")) + 
    geom_abline(intercept = 0, slope = 0, color = "blue", linewidth = .5) +
    geom_errorbar(aes(ymin=as.numeric(mean_dist)-as.numeric(sd_dist), ymax=as.numeric(mean_dist)+as.numeric(sd_dist))) +
    xlim(min(z_df$r2, 0), NA) + 
    ylim(min(z_df$mean_dist - as.numeric(max(as.numeric(z_df$sd_dist))), 0), max(z_df$mean_dist + as.numeric(z_df$sd_dist))) +
    theme_classic() +
    xlab(TeX("$\\r^2$"))
  
  return(p)
}

heatmap_gls <- function(z_df, phen, as = "as") {
  if (as == "as") {
    title <- paste0("Hairpin for ", phen, ", correcting for assessment center")
  } else {
    title <- paste0("Hairpin for ", phen, ", w/out correcting for assessment center")
  }
  
  ggplot(z_df, aes(factor(threshold), factor(pc_num), fill=on_line, width = 1, height = 1)) +
    geom_tile() +
    ggtitle(title) +
    theme_classic() +
    scale_fill_manual(values = c("on" = "black", "off" = "red")) +
    xlab("Threshold") +
    ylab("pc")
                         
}
