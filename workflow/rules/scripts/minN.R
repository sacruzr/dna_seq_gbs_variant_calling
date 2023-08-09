rm(list = ls())

pks <- c('dplyr','ggplot2','ggrepel','ggpubr','argparse')

pks <- suppressPackageStartupMessages(sapply(pks, require, character.only = TRUE))

if (any(!pks)) stop('The package(s) ', names(pks)[!pks], ' is/are not available for load.')

pr <- ArgumentParser()

pr$add_argument("-n", type = "integer", metavar = 'number',
                    help = "How many samples does the VCF have ?")

pr$add_argument("-i", type = "character", metavar = "file", default = "stdin",
                    help = "Input NS, one per line. No header")

pr$add_argument("-o", type = "character", metavar = "file",
                    help = "Output PDF file")

pr$add_argument("-t", type = "character", metavar = "title",
                    help = "Title of the plot", default = "NULL")

ar = pr$parse_args()

tNS <- ar$n
input <- ar$i
output <- ar$o
title <- ar$t

ns <-
  read.table(input, header = F) %>%
  rename(NS = V1) %>%
  group_by(NS) %>% 
  summarise(fNS = n()) %>% 
  full_join(x = data.frame(NS = 0:tNS), y = ., by = 'NS') %>% 
  replace(is.na(.), 0) %>% 
  mutate(c_vars = sum(fNS) - cumsum(fNS),
         c_calls = sum(NS * fNS) - cumsum(NS * fNS),
         p_miss = 1 - (c_calls / (c_vars * tNS))) %>% 
  replace(is.na(.), 0) %>%
  arrange(p_miss)


lab <- findInterval(seq(.1,1,.1), ns$p_miss)
names(lab) <-  seq(.1, 1, .1)
lab <-  lab[! (duplicated(lab) | duplicated(lab, fromLast = T))]

ns[lab,'labC'] <- paste( ns[lab,'NS'], 
                         scales::percent(ns[lab,'p_miss']),
                         sep = ' < NS\n' )

ns[lab,'labP'] <- paste( ns[lab,'NS'], 
                         format(ns[lab,'c_vars'], big.mark = ','),
                         sep = ' < NS\n' )

addUnits <- 
  function(n) {
  labels <- 
    ifelse(n < 1000, n,  # less than thousands
      ifelse(n < 1e6, paste0(round(n/1e3), 'k'),  # in thousands
        ifelse(n < 1e9, paste0(round(n/1e6), 'M'),  # in millions
          ifelse(n < 1e12, paste0(round(n/1e9), 'B'), # in billions
          'too big!'))))
  return(labels)
}

nx <- ceiling(max(ns$NS) * 0.075)
ny <- max(ns$p_miss) * 0.075

a <- 
  ggplot(data = ns, aes(x = NS, y = p_miss)) + 
  geom_point(size = 1, color = "#0F9D58") + theme_bw(base_size = 9) +
  geom_point(data = subset(ns, ! is.na(labP)), color = '#F4B400', size = 1.5) +
  labs(x = "Number of accessions genotyped (NS)", y = "Percentage of missing data") +
  geom_text_repel(aes(x = NS, y = p_miss, label = labP), nudge_x = nx, nudge_y = ny, size = 2) +
  scale_y_continuous(labels = scales::percent)

ny <- max(ns$c_vars) * 0.075

b <- 
  ggplot(data = ns, aes(x = NS, y = c_vars)) + 
  geom_point(size = 1, color = "#4285F4") + theme_bw(base_size = 9) +
  geom_point(data = subset(ns, ! is.na(labC)), color = '#F4B400', size = 1.5) +
  labs(x = "Number of accessions genotyped (NS)", y = "Number of variants") +
  geom_text_repel(aes(x = NS, y = c_vars, label = labC), nudge_x = nx, nudge_y = ny, size = 2) +
  scale_y_continuous(labels = addUnits)

c <- ggarrange(a, b, ncol = 2, nrow = 1)

if(title != 'NULL'){

	c <- annotate_figure(c, top = text_grob(title, face = "bold"))
}

ggsave(filename = output, plot = c, width = 170, height = 75, units = "mm")

