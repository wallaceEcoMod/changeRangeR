

#GET MASK FILES
spMaskDir=system.file("extdata/treeDemo/Migrationmask/120km/",package='changeRangeR')
allSpeciesMasks=tibble(rasterFiles=list.files(spMaskDir, recursive=T, full.names=T)) %>%
  mutate(sp.names= rasterFiles %>% basename %>% file_path_sans_ext)

maskCBSscenario(sumDirs = sumDirs,
                sp.ind = sp.ind,
                allSpeciesMasks = allSpeciesMasks,
                scenario = 'present',
                maskName = 'dispTry',
                mc.cores = 6)
