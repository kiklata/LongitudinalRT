# R O/E calculation
meta = seu@meta.data

meta.sub <- bb.res %>% dplyr::filter(project == "liver")

cluster.table <- table(meta.sub$celltype_major, meta.sub %>% rownames() %>% stringr::str_sub(.,18,18))

cluster.table %>%
  as.data.frame() %>%
  dplyr::mutate(
    p.value = purrr::pmap_dbl(
      list(
        .x = Var1,
        .y = Var2,
        .z = Freq
      ),
      .f = function(.x, .y, .z){
        a <- .z
        b <- sum(cluster.table[,.y]) - a
        c <- sum(cluster.table[.x,]) - a
        d <- sum(cluster.table) - a - b - c
        
        #o <- fisher.test(matrix(c(a, b,c, d), ncol = 2), alternative = "greater")
        #o$estimate
        o <- chisq.test(matrix(c(a, b, c, d), ncol = 2))
        oe <- o$observed/o$expected
        oe[1,1]
      }
    )
  ) -> enrich.res

#adj.p.value <- p.adjust(enrich.res$p.value, method = "BH")
#enrich.res <- enrich.res %>% dplyr::mutate(adj.p.value = adj.p.value)
enrich.res %>%
  dplyr::rename(`Ro/e` = p.value) %>%
  #dplyr::mutate(p.value = ifelse(p.value < 1, -1/p.value, p.value)) %>%
  #dplyr::mutate(`-log10(adj.P-value)` = -log10(adj.p.value)) %>%
  ggplot(aes(Var2, Var1, fill = `Ro/e`)) +
  geom_tile(colour = "white", lwd = 0.8) +
  theme(axis.title = element_text(size = 12)) +
  theme(axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 10)) +
  theme(legend.text = element_text(size = 10)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) +
  theme(axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color="black" ,angle = 45, hjust = 1)) +
  scale_fill_distiller(palette = "Spectral") -> roe.plot

roe.plot +  
  geom_text(data =  roe.plot$data %>% dplyr::filter(.,`Ro/e` > 1 & `Ro/e` < 2), 
            aes(label='*'), hjust = 0.5, vjust = 0.5)+
  geom_text(data =  roe.plot$data %>% dplyr::filter(.,`Ro/e` > 2 & `Ro/e` < 5), 
            aes(label='**'), hjust = 0.5, vjust = 0.5)+
  coord_equal()
                      