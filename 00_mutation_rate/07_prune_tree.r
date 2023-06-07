require(ape)

x <- read.nexus(file="oliveros_et_al_tree.txt")

tips_to_keep <- c("Tyrannus_albogularis", "Neopipo_cinnamomea", "Onychorhynchus_coronatus", "Pipra_filicauda")

x2 <- keep.tip(x, tips_to_keep)

plot(x2)
write.tree(x2, file="pruned_tree.tre")

x2 <- unroot(x2)

write.tree(x2, file="unrooted_tree.tre")
