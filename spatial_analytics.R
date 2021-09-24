#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

pkgs <- list('shiny' = '1.3.2',
             'shinyFiles' = '0.7.3',
             'data.table' = '1.12.2',
             'dplyr' = '0.8.3',
             'deldir' = '0.1-21',
             'doParallel' = '1.0.14',
             'foreach' = '1.4.4',
             'gplots' = '3.0.1.1',
             'ggplot2' = '3.2.1',
             'RColorBrewer' = '1.1-2',
             'tidyverse' = '1.2.1',
             'ComplexHeatmap' = '2.0.0',
             'circlize' = '0.4.6',
             'gtools' = '3.8.1',
             'spatstat' = '1.60-1',
             'Rfast' = '2.0.1',
             'ggpmisc' = '0.3.1')

msg <- function(r, p, v, o){
    x <- paste0('Please ', o, ' ', p, ' to at least version ', v, '.')
    if(is.na(r)){
        r <- x
    } else {
        r <- c(r, x)
    }
    return(r)
}

exit <- function() {
    .Internal(.invokeRestart(list(NULL, NULL), NULL))
}

instd <- as.data.frame(installed.packages(), stringsAsFactors = F)
report <- NA_character_

for( p in attributes(pkgs)$names){
    #print(p)
    if(p %in% instd$Package){
        if(packageVersion(p) < pkgs[[p]]){
            report <- msg(report, p, pkgs[[p]], 'update')
        } else {
            require(p, quietly = TRUE, character.only = TRUE)
        }
    } else {
        report <- msg(report, p, pkgs[[p]], 'install')
    }
}

if(!is.na(report)){
    print(paste(report, collapse = '\n'))
    exit()
}

options(shiny.maxRequestSize=4000*1024^2)

# Define UI for application that draws a histogram
ui <- fluidPage(
    # Application title
    titlePanel("Niche Finder"),

    sidebarPanel(
        # select output directory
        shinyDirButton("dir", "Chose output directory", "Upload"),

        textInput('anncol',
                  'Cluster Column Name',
                  'ClusterName'),

        textInput('efactcol',
                  'Patient / Treatment Group Column Name',
                  'Groups'),

        textInput('regcol',
                  'Column Name for Regions',
                  'Spots'),

        textInput('patcol',
                  'Column Name for Patients',
                  'Patients'),

        textInput('xcol',
                  'X Coordinate Column',
                  'X'),

        textInput('ycol',
                  'Y Coordinate Column',
                  'Y'),

        numericInput('numclust',
                     'Number of niches',
                     10,
                     min = 10,
                     max = 200),


        numericInput('mincells',
                     'Minimum cell count (total)',
                     100,
                     min = 50,
                     max = 500),

        numericInput('convfactor',
                     'Conversion factor (um/pixel)',
                     0.37744,
                     min = 0.05,
                     max = 5),
        radioButtons('denshis', 'Distance graphs:',
                     c('density' = 'dens',
                     'histogram' = 'hist'),
                     inline = TRUE),

        # Input: Select a file ----
        fileInput(
            "data_file",
            "Choose CSV File",
            multiple = FALSE,
            accept = c("text/csv",
                       "text/comma-separated-values,text/plain",
                       ".csv")
        )
    ),

    mainPanel(

        verbatimTextOutput('path'),
        br(),
        textOutput('delmsg'),
        textOutput('pivotmsg'),
        textOutput('llhmsg'),
        textOutput('eqalmsg'),
        textOutput('globalmsg'),
        textOutput('nichemsg'),
        textOutput('profmsg'),
        textOutput('celldistmsg')
    )
)


# Define server logic required to draw a histogram
server <- function(input, output) {

    x_coord_colname <- 'X:X'
    y_coord_colname <- 'Y:Y'
    ann_colname     <- 'celltypes'
    tile_colname    <- 'region'
    efact_colname   <- 'efact' # experimental factor colname
    #home_dir        <- '~/mnt/driveB/jdemeter/project/nolanlab/'
    #home_dir        <- 'E:/CHRISTIAN/For Janos/sizun/'
    home_dir        <- 'z:/project/nolanlab/'
    #cat('homedir = ', home_dir)
    perm_no         <- 1000
    folders         <- list(del = 'delauney_files',
                            int = 'intermediate_files',
                            niche = 'niche_heatmaps',
                            ct = 'celltype_heatmaps',
                            dist = 'celldist_plots')

    concat <- function(a, b, sep = '__'){
        a <- unlist(a)
        b <- unlist(b)
        return(paste0(a, sep, b))
    }

    clean_path <- function(txt){
        gsub("[^a-zA-Z0-9\\.\\-]", "_", txt)
    }

    make_matr <- function(ay, sep = '___'){
        matrix(
            apply(expand.grid(ay, ay),
                  1,
                  paste, collapse = sep),
            nrow = length(ay)
        )
    }

    mod_dset <- function(inp){

        dt <- fread(inp$data_file$datapath) %>%
            mutate(!!x_coord_colname := get(inp$xcol),
                   !!y_coord_colname := get(inp$ycol),
                   !!tile_colname := as.character(get(inp$regcol)),
                   !!ann_colname := get(inp$anncol)) %>%
            as.data.table

        if(inp$efactcol != 'none' & inp$efactcol %in% names(dt)){
            dt <- dt %>%
                mutate(!!efact_colname := get(inp$efactcol)) %>%
                as.data.table
        } else {
            dt <- dt %>%
                mutate(!!efact_colname := 1) %>%
                as.data.table
        }

        if (!'XYcellid' %in% colnames(dt)) {
            if (all(c(x_coord_colname, y_coord_colname) %in% colnames(dt))){
                dt[, XYcellid := concat(dt[[x_coord_colname]], dt[[y_coord_colname]])]
            }
        }

        dt[!grepl('^reg', get(tile_colname)), (tile_colname) := paste0('reg', get(tile_colname))]

        dt <- dt %>%
            select(tile_colname, x_coord_colname, y_coord_colname, ann_colname,
                   XYcellid, efact_colname) %>%
            as.data.table

        # filter out regions with only 1 cell
        if(nrow(dt[, .N, tile_colname][N==1]) > 0){
            singlets <- dt[, .N, tile_colname][N==1, ][[tile_colname]]
            dt <- dt[!get(tile_colname) %in% singlets]
            cat(file = stderr(), paste('removed tile(s):', singlets, '(singlets)\n'))
        }
        # filter out regions with 2 cells that share the same x/y coordinate
        if(nrow(dt[, .N, tile_colname][N==2]) > 0){
            bad_tiles <- dt[, .N, tile_colname][N == 2, ][[tile_colname]] %>% unlist %>% unname
            for(t in bad_tiles){
                if(length(unique(dt[get(tile_colname) == t, `X:X`])) == 1 |
                   length(unique(dt[get(tile_colname) == t, `Y:Y`])) == 1){
                    dt <- dt[get(tile_colname) != t]
                    cat(file = stderr(), paste('removed tile:', t,'(clashing coords for 2 cells)\n'))
                }
            }
        }

        # if regions are names like regNNN, make sure we don't have leading 0s
        if(all(grepl('^reg\\d+$', dt[[tile_colname]]))){
            dt[, (tile_colname) := lapply(.SD,
                                          function(x){
                                              paste0('reg', as.integer(gsub('^reg(\\d+)$', '\\1', x)))
                                          }),
               .SDcols = tile_colname]
        }

        # make sure each cell is unique:
        dt <- dt[, .SD[1,], .(region, XYcellid)]

        #cat(file=stderr(), paste('before: nrow=', nrow(dt), ',  celltypes=', unique(dt$celltypes), '\n'))
        # filter out cell types, that have cell numbers < mincells
        dt <- dt[celltypes %in% dt[, .N, c('celltypes', efact_colname)][N > min_cells(), celltypes]]
        #cat(file=stderr(), paste('after: nrow=', nrow(dt), ',  celltypes=', unique(dt$celltypes), '\n'))

        cat(file=stderr(), 'Imported dataset ...\n')
        cat(file=stderr(), file.path(path(), "dt.txt\n"))

        # cleanup folder
        # unlink(paste(path(), '*', sep = '/'), recursive = TRUE)

        # create subfolders
        for(f in folders){
            dir.create(file.path(path(), f))
        }

        fwrite(dt, file.path(path(), "dt.txt"), sep = '\t')
        return(dt)
    }

    # dir
    # cat('homedir=', home_dir)
    shinyDirChoose(input, 'dir', roots = c(home = home_dir), filetypes = c('', 'txt'))
    out_dir    <- reactive({input$dir})

    # path
    path <- reactiveVal(getwd())
    observeEvent(ignoreNULL = TRUE,
                 eventExpr = {
                     input$dir
                 },
                 handlerExpr = {
                     if (!"path" %in% names(out_dir())) return()
                     home <- normalizePath(home_dir)
                     path(file.path(home, paste(unlist(out_dir()$path[-1]), collapse = .Platform$file.sep)))
                     # create dir inside the selected folder
                     path(file.path(path(),
                                    paste0('nf_',
                                           gsub('[^0-9]', '', format(Sys.time(), "%d%m%y%X")))))
                     dir.create(path())
                 })

    output$dir  <- renderPrint(path)

    annot_col   <- reactive({input$anncol})
    tissue_col  <- reactive({input$regcol})
    efact_col   <- reactive({input$efactcol})
    x_col       <- reactive({input$xcol})
    y_col       <- reactive({input$ycol})
    pat_col     <- reactive({input$patcol})
    inputf_ori  <- reactive({input$data_file$datapath})
    dt          <- reactive({mod_dset(input)})
    output$path <- renderPrint(path())
    no_clust    <- reactive({input$numclust})
    min_cells   <- reactive({input$mincells})
    conv_factor <- reactive({input$convfactor})
    dens_his    <- reactive({input$denshis})

    output$delmsg <- renderPrint({

        req(input$data_file)
        dt        <- dt()
        outDir    <- path()
        annotCol  <- ann_colname  #annot_col()
        tissueCol <- tile_colname  #tissuve_col()
        tissues   <- dt[[tissueCol]] %>% unique()

        cores <- detectCores()
        cl    <- makeCluster(min(c(length(tissues), cores)),
                             outfile = file.path(outDir, folders$int, 'dopar_log'))
        registerDoParallel(cl)

        summary_file <- file.path(outDir, folders$del, "deldir_summary.tsv")
        empty_dt <- data.table(x = integer(), y = integer(), n.tri = integer(), del.area = numeric(),
                               del.wts = numeric(), n.tside = integer(), nbpt = integer(),
                               dir.area = numeric(), dir.wts = numeric())
        #file.create(summary_file)
        fwrite(empty_dt, summary_file, sep = '\t')
        re <<- foreach(
            i = 1:length(tissues),
            .combine = 'bind_rows',
            .export = c('concat'),
            .packages = c("deldir", 'data.table', 'dplyr')
        ) %dopar% {
            print(i)
            sub    <- dt[get(tissueCol) == tissues[i],]
            tryCatch({
                vtess  <- deldir(sub[, `X:X`], sub[, `Y:Y`])
                fwrite(vtess$summary %>% mutate(tissueCol = tissues[i]) %>% as.data.table,
                       summary_file, append = TRUE, sep = '\t')
                result <- vtess$delsgs %>% select(-ind1,-ind2) %>% as.data.table
                result[, cell1ID := concat(result[, 1], result[, 2])]
                result[, cell2ID := concat(result[, 3], result[, 4])]

                result <- result %>% left_join(sub %>% select(annotCol, XYcellid),
                                               by = c('cell1ID' = 'XYcellid')) %>% as.data.table
                result <- result %>% left_join(sub %>% select(annotCol, XYcellid),
                                               by = c('cell2ID' = 'XYcellid')) %>% as.data.table

                result <- result[, c(1:4, 7:8)]
                colnames(result)[c(5, 6)] = c("cell1type", "cell2type")
                result <- rbind(result, result[, .(x1 = x2, y1 = y2, x2 = x1, y2 = y1,
                                                   cell1type = cell2type, cell2type = cell1type)])
                print('done')
                result[, ind := i]
            }, error = function(err){
                data.table(x1 = integer(), y1 = integer(), x2 = integer(), y2 = integer(),
                           cell1type = character(), cell2type = character(), ind = integer())
            })
        }
        stopCluster(cl)
        for (i in seq(length(tissues))) {
            fname <- file.path(outDir, folders$del, paste(tissues[i], 'Rdelaun.csv', sep = '_'))
            fwrite(re[ind == i, -'ind'], fname)
        }
        cat(file=stderr(), 'Finished Delauney graph calculations ... \n')
        print('Finished Delauney graph calculations')
    })

    output$pivotmsg <- renderPrint({

        req(input$data_file)
        outDir <- path()

        for (fname in dir(file.path(outDir, folders$del), full.names = TRUE)) {
            if (grepl('^reg.+Rdelaun.csv', basename(fname))) {
                #print(fname)
                delaun <- fread(fname, header = TRUE)

                pivot <- delaun %>% dcast(cell1type ~ cell2type,
                                          value.var = 'x1',
                                          fun.aggregate = length) %>%
                    as.data.table

                # normalize values by dividing with row sums
                x <- rowSums(pivot[, 2:ncol(pivot)])
                pivotnorm <- bind_cols(pivot[, .(cell1type)],
                                       sapply(2:ncol(pivot),
                                              function(i) {
                                                  pivot[[i]] / x[i - 1]
                                              }) %>% as.data.table)

                orifname  <- file.path(outDir, folders$int, paste0("pivotORI_", basename(fname)))
                normfname <- file.path(outDir, folders$int, paste0("pivotNORM_", basename(fname)))
                fwrite(pivot, orifname)
                fwrite(pivotnorm, normfname)
            }
        }
        cat(file=stderr(), 'Done with pivot tables ...\n')
        print('Finished creating pivot tables.')
    })

    output$llhmsg <- renderPrint({

        req(input$data_file)
        outDir <- path()

        for (fname in dir(file.path(outDir, folders$int), full.names = TRUE)) {
            if (grepl('^pivotORI_reg.*Rdelaun.csv$', basename(fname))) {

                #reading input file
                pivot <- fread(fname)

                # remove dirt column/row
                pivot <- pivot[cell1type != "dirt", ]
                if('dirt' %in% colnames(pivot)){
                    pivot <- pivot[, -'dirt']
                }

                # make matrix from data
                ma    <- as.matrix(pivot[, -'cell1type'])
                rownames(ma) <- pivot$cell1type
                colnames(ma) <- colnames(pivot)[-1]

                # normalize values by rowsum and totalsums
                summa   <- sum(ma)
                hornorm <- unname(rowSums(ma)) / summa

                for (i in seq(nrow(ma))) {
                    for (j in seq(ncol(ma))) {
                        ma[i, j] = ma[i, j] / (hornorm[i] * hornorm[j] * summa)
                    }
                }

                ma <- as.data.table(ma, keep.rownames = 'pairs')

                llhfname <- file.path(outDir,
                                      folders$int,
                                      paste("lglikelihood_mx_from", basename(fname), sep = '_'))
                fwrite(ma, llhfname)
            }
        }
        cat(file=stderr(), 'Finished likelihood calculations ...\n')
        print('finished lglikelihood calculations.')
    })

    output$eqalmsg <- renderPrint({

        req(input$data_file)
        cluster_or_reorder <- 'reorder'
        outDir             <- path()
        combo              <- data.table()

        for(fname in dir(file.path(outDir, folders$int), full.names = TRUE)){
            if(grepl('^lglikelihood_mx_from_pivotORI_reg.+_Rdelaun.csv$', basename(fname))){

                fstub  <- gsub('\\.csv', '', basename(fname))
                llhtab <- fread(fname)
                a <- as.matrix(llhtab[, -1])
                ##################### done reading and re-formatting the matrix
                b <- make_matr(colnames(llhtab)[-1])
                #####################done making the matrix of row and column names for the interaction matrix
                b <- b[!lower.tri(b)]
                ##################done unfolding the row names for the upper triangle of the interaction matrix
                a <- a[!lower.tri(a)]
                ##################done unfolding the values for the upper triangle of the interaction matrix
                if (nrow(combo) == 0){
                    combo <- data.table(annotation = b, fname = a)
                    colnames(combo)[2] <- fstub
                    ################done naming a row after file
                } else {
                    c <- data.table(annotation = b, fname = a)
                    colnames(c)[2] <- fstub
                    combo <- full_join(combo, c, by = "annotation") %>% as.data.table
                }
            }
        }

        combo[, leftcell := gsub('^(.+)___.+$', '\\1', annotation)]
        combo[, rightcell := gsub('^(.+)___(.+)$', '\\2', annotation)]
        setcolorder(combo, c('leftcell', 'rightcell', colnames(combo)[!colnames(combo) %in% c('leftcell', 'rightcell')]))

        combo     <- unique(combo)

        breaks    <- seq(0, 4, length.out = 1000)
        gradient1 <- colorpanel( sum( breaks[-1]<=1 ), "#0030FF", "white" )
        gradient2 <- colorpanel( sum( breaks[-1]>1 ), "white", "#FF5020" )
        hm.colors <- c(gradient1,gradient2)

        dd        <- function(c) as.dist(1 / (c + 5)) # from Yury
        data_ref  <- function(dt, colind, rpl = F){
            col <- colnames(dt)[colind]
            z <- dt %>%
                dcast(leftcell ~ rightcell, value.var = col, fill = 1) %>%
                as.data.table %>%
                .[order(leftcell)] %>%
                setcolorder(.[['leftcell']]) %>%
                as.data.frame %>%
                column_to_rownames('leftcell') %>%
                as.matrix()

            if(rpl == TRUE){
                z[is.na(z)] <- 0
            }
            return( z )
        }

        if (cluster_or_reorder =="reorder"){
            print ("reordering")

            # order of rows is based on one of the files (one with no
            # missing values)
            colind   <- min(
                which(combo[, lapply(.SD[, -c(1:3)],
                                     function(x){
                                         sum(is.na(x))
                                     })] == 0)) + 3
            z        <- data_ref(combo, colind, rpl = T)
            # safer to print to disposable file
            png(filename = file.path(outDir, folders$ct, 'waste.png'),
                height = 1500, width = 2000, res = 300, pointsize = 10)
            hmstore  <- heatmap.2(z, distfun = dd, Rowv = F, dendrogram = 'column')
            dev.off()
            unlink(file.path(outDir, folders$ct, 'waste.png'))
        }

        for(i in seq(4, ncol(combo))){

            z <- data_ref(combo, i)

            pngfile     <- file.path(outDir, folders$ct, paste(colnames(combo)[i], 'png', sep='.'))
            png(filename = pngfile,height = 1500, width = 2000, res = 300, pointsize = 10)
            if(cluster_or_reorder =="reorder"){
                #z       <- z[hmstore$rowInd, hmstore$rowInd]
                heatmap.2(z, col=hm.colors, trace="none", lmat=rbind(4:3,2:1),
                          breaks=breaks, margins=c(15,15), cexRow=0.5,
                          cexCol=0.5, lwid=c(1.5,2.0), Colv = FALSE, Rowv = FALSE,
                          dendrogram = 'none', na.color = 'grey60')
            } else {
                heatmap.2(z, distfun=dd, col=hm.colors, trace="none",
                          lmat=rbind(4:3,2:1), breaks=breaks, margins=c(15,15),
                          cexRow=0.5, cexCol=0.5, lwid=c(1.5,2.0),
                          na.color = 'grey60')
            }
            dev.off()
        }
        cat(file = stderr(), 'Ready with clustering of llh by region ...\n')
        print('ready with clustering of llh by region')
    })

    output$globalmsg <- renderPrint({

        # create heatmap of cluster1 vs cluster2 using likelihood ratio and frequency metrics
        # also, do this by groups (samples = values of expt factor)

        req(input$data_file)

        outDir <- path()
        dtmap  <- dt() %>%
            select(tile_colname, efact_colname) %>%
            unique %>%
            as.data.table

        combo  <- data.table()

        for(fname in dir(file.path(outDir, folders$del), full.names = TRUE)){
            if(grepl('^reg.+Rdelaun.csv', basename(fname))){
                d <- fread(fname)
                region <- gsub('^(reg\\d+)_.*$', '\\1', basename(fname))
                d[, reg := region]
                combo <- bind_rows(combo, d)
            }
        }

        combo <- combo %>%
            left_join(dtmap, by = c('reg' = 'region')) %>%
            as.data.table

        comb_res <- function(a, b){
            full_join(a, b, by = c('cell1type', 'cell2type', 'efact')) %>% as.data.table
        }

        scramble_cols <- function(dt, efs, n = perm_no){

            cores <- detectCores()
            cl    <- makeCluster(cores, outfile = file.path(outDir, folders$int, 'dopar_log'))
            registerDoParallel(cl)

            set.seed <- 12345

            re <- foreach(
                i = 1:n,
                .combine = 'comb_res',
                .inorder = FALSE,
                .export = c('calc_metrics'),
                .packages = c('data.table', 'dplyr')
            ) %dopar% {
                dtx <- copy(dt)
                for(ef in efs){
                    dtx[efact == ef, `:=`(cell1type = sample(cell1type, replace = T), cell2type = sample(cell2type, replace = T))]
                }
                dtx <- calc_metrics(dtx, flag = 2)
                dtx <- dtx[, .(cell1type, cell2type, efact, llh)]
                #print(dtx)
                dtx
            }
            stopCluster(cl)
            colnames(re)[4:ncol(re)] <- paste0('llh', seq(1, ncol(re)-3))

            # make sure we have all permutations or celltypes for all efacts
            full_re <- data.table()
            for(ef in efs){
                full_re <- bind_rows(full_re, unique(dt$cell1type) %>% permutations(length(.), 2, ., repeats.allowed = T) %>%
                    as.data.table %>%
                    unite('key', V1, V2, sep = '__', remove = FALSE) %>%
                    rename('cell1type' = V1, 'cell2type' = V2) %>%
                    mutate(efact = ef) %>%
                    as.data.table)
            }

            re <- full_re %>%
                full_join(re %>% unite('key', cell1type, cell2type, sep = '__', remove = FALSE),
                          by = c('key', 'cell1type', 'cell2type', 'efact')) %>%
                select(-key) %>%
                as.data.table

            invisible(re)
        }

        calc_metrics <- function(dt, flag = 1){
            dtx <- dt[, .(ints_cl12 = .N), .(cell1type, cell2type, efact)] %>%
                left_join(dt[, .(ints_cl1 = .N), .(cell1type, efact)],
                          by = c('cell1type', 'efact')) %>% # interaction counts per cell1
                left_join(dt[, .(ints_cl2 = .N), .(cell2type, efact)],
                          by = c('cell2type', 'efact')) %>% # interaction counts per cell2
                left_join(dt[, .(ints_total = .N), efact], by = 'efact') %>% # all interactions
                # calculate likelihood ratio for interactions
                mutate(llh = (1.0 * ints_cl12 * ints_total) / (1.0 * ints_cl1 * ints_cl2)) %>%
                as.data.table
            if(flag == 1){
                dtx <- dtx %>%
                    # count cells
                    left_join(dt[, .(cell_cl12 = nrow(unique(.SD[, .(x1,y1,region)]))), .(cell1type, cell2type, efact)],
                              by = c('cell1type', 'cell2type', 'efact')) %>%
                    left_join(dt[, .(cell_cl21 = nrow(unique(.SD[, .(x2,y2,region)]))), .(cell1type, cell2type, efact)],
                              by = c('cell1type', 'cell2type', 'efact')) %>%
                    left_join(dt[, .(cell_cl1 = nrow(unique(.SD[, .(x1,y1,region)]))), .(cell1type, efact)],
                              by = c('cell1type', 'efact')) %>%
                    left_join(dt[, .(cell_cl2 = nrow(unique(.SD[, .(x2,y2,region)]))), .(cell2type, efact)],
                              by = c('cell2type', 'efact')) %>%
                    # calculate interaction frequency
                    mutate(`cl12_cl1` = 1.0 * ints_cl12/ints_cl1) %>%
                    as.data.table
            }
            invisible(dtx)
        }

        combox <- calc_metrics(combo)
        fwrite(combox, file.path(outDir, folders$ct, 'cell1type_vs_cell2type_llh.txt'), sep = '\t')

        efacts <- sort(unique(combox$efact))

        llh_bg <- scramble_cols(combo, efacts)
        fwrite(llh_bg, file.path(outDir, folders$ct, 'scrambled_data.txt'), sep = '\t')

        cell_count_max <- round(max(combox$cell_cl1))
        llh_max <- max(abs(range(log2(combox$llh)))) * 0.7
        llh_min <- -1 * llh_max
        mats <- list()
        dfs  <- list()

        for (ef in efacts) {

            # calculate p-values
            z <- combox[efact == ef, .(cell1type, cell2type, llh)] %>%
                left_join(llh_bg[efact == ef, -'efact'], by = c('cell1type', 'cell2type')) %>%
                as.data.table

            z[, mean := apply(.SD, 1, mean, na.rm = T), .SDcols = grep('^llh\\d', colnames(z))]
            z[, sd := apply(.SD, 1, sd, na.rm = T), .SDcols = grep('^llh\\d', colnames(z))]
            z[, qt := pnorm(llh, mean, sd)]
            z[, pval1 := pmin(2-2*qt, 2*qt, na.rm = T)]
            z[, padj1 := p.adjust(pval1, method = 'bonferroni')]
            z[, pval := apply(.SD, 1, function(x){
                if(sum(!is.na(x[-1]))>0){
                    F <- ecdf(x[-1])
                    min(F(x[1])*2, 2-2*F(x[1]))
                } else {
                    NA_real_
                }}), .SDcols = grep('^llh', colnames(z))]
            z[, padj := lapply(.SD, p.adjust, method = 'BH'), .SDcols = 'pval']
            # z <- z[, .(cell1type, cell2type, padj)]
            # pvmin <- -1*(z$padj[z$padj > 0]%>% log10 %>% min %>% floor)

            # pvalues in this matrix
            zmat <- z %>%
                dcast(cell1type ~ cell2type, value.var = 'padj1') %>%
                melt(id.var = 'cell1type',
                     variable.name = 'cell2type',
                     value.name = 'padj') %>%
                mutate(lpadj = case_when(!is.na(padj) ~ -log10(padj),
                                         #padj > 0 ~ -log10(padj),
                                         #padj == 0 ~ -log10(padj), #pvmin,
                                         is.na(padj) ~ 0.)) %>%
                as.data.table %>%
                #mutate(tpadj = as.character(round(lpadj,1))) %>%
                dcast(cell1type ~ cell2type, value.var = 'lpadj') %>%
                dplyr::select(sort(colnames(.))) %>%
                arrange(cell1type) %>%
                as.data.frame %>%
                column_to_rownames('cell1type') %>%
                #select(noquote(colnames(.))) %>% # This now gives an error
                as.matrix

            # p-values as text in this matrix
            zmatt <- z %>%
                dcast(cell1type ~ cell2type, value.var = 'padj1') %>%
                melt(id.var = 'cell1type',
                     variable.name = 'cell2type',
                     value.name = 'padj') %>%
                mutate(lpadj = case_when(!is.na(padj) ~ -log10(padj),
                                         #padj > 0 ~ -log10(padj),
                                         #padj == 0 ~ -log10(padj), #pvmin,
                                         is.na(padj) ~ 0.)) %>%
                mutate(tpadj = as.character(round(lpadj, 2))) %>%
                as.data.table %>%
                dcast(cell1type ~ cell2type, value.var = 'tpadj') %>%
                dplyr::select(sort(colnames(.))) %>%
                arrange(cell1type) %>%
                as.data.frame %>%
                column_to_rownames('cell1type') %>%
                #select(noquote(colnames(.))) %>% # this now gives an error
                as.matrix

            # likelihood ratio values
            mat <- combox[efact == ef,] %>%
                dcast(cell1type ~ cell2type, value.var = 'llh') %>%
                melt(id.var = 'cell1type',
                     variable.name = 'cell2type',
                     value.name = 'llh') %>%
                mutate(lllh = case_when(llh > 0 ~ log2(llh),
                                        is.na(llh) ~ 0)) %>%
                as.data.table %>%
                dcast(cell1type ~ cell2type, value.var = 'lllh') %>%
                dplyr::select(sort(colnames(.))) %>%
                arrange(cell1type) %>%
                as.data.frame %>%
                column_to_rownames('cell1type') %>%
                #select(noquote(colnames(.))) %>% # this now gives an error
                as.matrix

            mats[[length(mats)+1]] <- mat

            # cell counts for top bar
            df <- unique(combox[efact == ef, .(cell_count = cell_cl1, cell1type)]) %>%
                as.data.frame(stringsAsFactors = FALSE) %>%
                right_join(data.frame(cell1type = rownames(mat), stringsAsFactors = FALSE), by = 'cell1type') %>%
                arrange(cell1type) %>%
                column_to_rownames('cell1type')

            # top color bar for cell numbers
            ha <- HeatmapAnnotation(df = df,
                                    col = list(cell_count = colorRamp2(c(0, cell_count_max),
                                                                       c("white", "red"))))

            dfs[[length(dfs)+1]] <- df

            ef_postfix <- ifelse(length(efacts)==1 & ef==1, '', paste0('_', clean_path(ef)))
            pdffile <- file.path(outDir, folders$ct, paste0('cell1type_vs_cell2type_llh', ef_postfix, '.pdf'))
            pdf(file = pdffile, width = 11, height = 10)
            print(
                mat %>%
                    Heatmap(
                        cluster_columns = FALSE,
                        cluster_rows = FALSE,
                        na_col = 'grey90',
                        col = circlize::colorRamp2(seq(from = llh_min, to = llh_max, length.out = 100),
                                                   colorRampPalette(rev(brewer.pal(n= 7, 'RdYlBu')))(100)),
                        top_annotation = ha,
                        column_title = paste('likelihood ratios for sample:', ef),
                        heatmap_legend_param = list(title = 'log2 lhr', at = c(round(llh_min), 0, round(llh_max))),
                        # add pval to graph
                        cell_fun = function(j, i, x, y, width, height, fill) {
                            if(!is.na(zmat[i, j]) & zmat[i,j] > -log10(0.01))
                                grid.text(zmatt[i, j], x, y, gp = gpar(fontsize = 5))
                        }
                    )
            )
            dev.off()
        }

        ##### plot comparison plots
        matdiffs <- list()
        matps    <- list()
        matpst   <- list()
        dfdiffs  <- list()
        set.seed(12345)
        if(length(efacts) > 1){
            comparisons <- combinations(length(efacts), 2, v = efacts)
            for(r in 1:nrow(comparisons)){
                a <- comparisons[r,1]
                b <- comparisons[r,2]
                mata <- mats[[which(efacts == a)]]
                matb <- mats[[which(efacts == b)]]

                isect <- intersect(rownames(mata), rownames(matb))
                mata  <- subset(mata, rownames(mata) %in% isect) %>% .[, isect]
                matb  <- subset(matb, rownames(matb) %in% isect) %>% .[, isect]

                matdiffs[[r]] <- matb - mata
                dfdiffs[[r]] <- dfs[[which(efacts == a)]] %>% as.data.table(keep.rownames = T) %>%
                    inner_join(dfs[[which(efacts == b)]] %>% as.data.table(keep.rownames = T),
                              by = 'rn',
                              suffix = c(paste0('_', a), paste0('_', b))) %>%
                    column_to_rownames('rn')

                # get distr of randomized log ratios
                difllh <- matb %>%
                    as.data.table(keep.rownames = T) %>%
                    melt(id.vars = 'rn') %>%
                    left_join(mata %>%
                                  as.data.table(keep.rownames = T) %>%
                                  melt(id.vars = 'rn'),
                              by = c('rn', 'variable'),
                              suffix = c('.b', '.a')) %>%
                    mutate(llh = value.b - value.a) %>%
                    select(-value.a, -value.b) %>%
                    unite('key', rn, variable, sep = '__') %>%
                    arrange(key) %>%
                    as.data.table

                difllh <- difllh %>% left_join(
                    (llh_bg[efact == b & cell1type %in% isect & cell2type %in% isect, -'efact'] %>%
                         unite('rn', cell1type, cell2type, sep = '__') %>%
                         arrange(rn) %>%
                         column_to_rownames('rn') %>%
                         as.matrix %>%
                         log2
                     -
                     llh_bg[efact == a & cell1type %in% isect & cell2type %in% isect, -'efact'] %>%
                         unite('rn', cell1type, cell2type, sep = '__') %>%
                         arrange(rn) %>%
                         column_to_rownames('rn') %>%
                         as.matrix %>%
                         log2
                     ) %>%
                        as.data.table(keep.rownames = TRUE), by = c('key' = 'rn')) %>%
                    as.data.table

                difllh[, mean := apply(.SD, 1, mean, na.rm = T), .SDcols = grep('^llh\\d', colnames(difllh))]
                difllh[, sd := apply(.SD, 1, sd, na.rm = T), .SDcols = grep('^llh\\d', colnames(difllh))]
                difllh[, qt := pnorm(llh, mean, sd)]
                difllh[, pval1 := pmin(2-2*qt, 2*qt, na.rm = T)]
                difllh[, padj1 := p.adjust(pval1, method = 'bonferroni')]
                difllh[, pval := apply(.SD, 1, function(x){
                    if(sum(!is.na(x[-1]))>0){
                        F <- ecdf(x[-1])
                        min(F(x[1])*2, 2-2*F(x[1]))}
                    else{
                        NA_real_
                    }}), .SDcols = grep('^llh', colnames(difllh))]
                difllh[, padj := lapply(.SD, p.adjust, method = 'BH'), .SDcols = 'pval']

                matps[[r]] <- difllh[, .(key, padj = padj1)] %>%
                    separate(key, into = c('cell1type', 'cell2type'), sep = '__') %>%
                    mutate(lpadj = case_when(is.na(padj) ~ 0.,
                                             !is.na(padj) ~ -log10(padj))) %>%
                    dcast(cell1type ~ cell2type, value.var = 'lpadj') %>%
                    as.data.table %>%
                    dplyr::select(sort(colnames(.))) %>%
                    arrange(cell1type) %>%
                    as.data.frame %>%
                    column_to_rownames('cell1type') %>%
                    as.matrix

                matpst[[r]] <- difllh[, .(key, padj = padj1)] %>%
                    separate(key, into = c('cell1type', 'cell2type'), sep = '__') %>%
                    mutate(lpadj = case_when(is.na(padj) ~ 0.00,
                                             !is.na(padj) ~ -log10(padj))) %>%
                    mutate(tpadj = as.character(round(lpadj, 2))) %>%
                    as.data.table %>%
                    dcast(cell1type ~ cell2type, value.var = 'tpadj') %>%
                    dplyr::select(sort(colnames(.))) %>%
                    arrange(cell1type) %>%
                    as.data.frame %>%
                    column_to_rownames('cell1type') %>%
                    as.matrix
                #
                # diffba <- ecdf(
                #     (sample(as.numeric(matb), size = 100000, replace = T)
                #      -
                #      sample(as.numeric(mata), size = 100000, replace = T)))
                #
                # adjusted pvals for observed log-ratios
                # matps[[r]] <- matrix(
                #     #p.adjust(
                #         sapply(as.numeric(matdiffs[[r]]),
                #                function(x){
                #                    min(c(diffba(x)*2, (1-diffba(x))*2))
                #                }),
                #     #    method = 'BH'),
                #     nrow = nrow(mata)) # %>% log10
                #
                #matdiffs[[r]]
                #colnames(dfdiffs[[length(dfdiffs)]]) <- c(a,b)
            }

            gmin <- min(sapply(matdiffs, min))
            gmax <- max(sapply(matdiffs, max))

            for(k in 1:nrow(comparisons)){

                colmap <- list()
                cname1 <- colnames(dfdiffs[[k]])[1]
                cname2 <- colnames(dfdiffs[[k]])[2]
                colmap[[cname1]] <- colorRamp2(c(0, cell_count_max), c("white", "red"))
                colmap[[cname2]] <- colorRamp2(c(0, cell_count_max), c("white", "red"))
                legmap <- list()
                legmap[[cname1]] <- list(title = 'Cell Count')
                legmap[[cname2]] <- list(title = 'Cell Count')

                ha <- HeatmapAnnotation(df = dfdiffs[[k]], col = colmap,
                                        annotation_legend_param = legmap,
                                        show_legend = c(TRUE, FALSE))

                pdffile <- file.path(outDir, folders$ct, paste0('cell1type_vs_cell2type_x',
                                                                paste(comparisons[k,2],
                                                                      comparisons[k,1],
                                                                      sep = '_vs_'), '.pdf'))
                pdf(file = pdffile, width = 11, height = 10)
                print(
                    matdiffs[[k]] %>%
                        Heatmap(
                            cluster_columns = TRUE,
                            cluster_rows = TRUE,
                            na_col = 'grey90',
                            col = circlize::colorRamp2(seq(from = gmin, to = gmax, length.out = 100),
                                                       colorRampPalette(rev(brewer.pal(n= 7, 'BrBG')))(100)),
                            top_annotation = ha,
                            column_title = paste('log2fc for likelihood ratios:',
                                                 comparisons[k,2], 'vs.', comparisons[k,1]),
                            heatmap_legend_param = list(title = 'l2fc lhr'),
                            rect_gp = gpar(col = "grey80", lwd = 1),
                            cell_fun = function(j, i, x, y, width, height, fill) {
                                if(!is.na(matps[[k]][i, j]) & matps[[k]][i,j] > -log10(0.001)){
                                    grid.text(matpst[[k]][i, j], x, y, gp = gpar(fontsize = 5))
                            }}
                        )
                )
                dev.off()
            }

        }
        cat(file=stderr(), 'Ready with clustering of llh globally ...\n')
        print('ready with clustering of llh globally')
    })

# niches --------------------------------------------------------------------------------------

    output$nichemsg <- renderPrint({

        req(input$data_file)
        outDir      <- path()
        dt          <- dt()
        annotCol    <- ann_colname  #annot_col()
        tissueCol   <- tile_colname  #tissue_col()
        efactCol    <- efact_colname

        ef2reg_map  <- dt[, .(region, efact)] %>% unique %>% .[order(efact, region)]
        efacts      <- unique(ef2reg_map$efact)

        upperTrivec <- function(dtab){
            dtab <- as.matrix(dtab)
            sapply(1:(nrow(dtab)), function(i){dtab[i,seq(i, nrow(dtab), 1)]}) %>% unlist
        }

        for(j in 1:length(efacts)){

            reg_counts <- dcast(dt[get(efactCol) == efacts[j]],
                                as.formula(paste(annotCol, '~', tissueCol)),
                                value.var = 'X:X', fun.aggregate = length)

            setcolorder(reg_counts,
                        c(annotCol,
                          paste0('reg', sort(as.integer(gsub('reg', '', colnames(reg_counts)[-1]))))))

            tissues <- dt[get(efactCol) == efacts[j],][[ tissueCol ]] %>%
                unlist %>%
                unique %>%
                str_replace('reg', '') %>%
                as.integer %>%
                sort() %>%
                paste0('reg', .)

            combo <- data.table()

            for (i in seq(1, length(tissues))){

                #print(tissues[i])

                #reading lglilyhd_input file
                lglklyhd_input <- fread(
                    file.path(
                        outDir,
                        folders$int,
                        paste("lglikelihood_mx_from_pivotORI", tissues[i], "Rdelaun.csv", sep = "_")
                )) %>%
                    filter(pairs != 'dirt') %>%
                    select(colnames(.)[colnames(.) != 'dirt']) %>%
                    as.data.frame %>%
                    column_to_rownames(colnames(.)[1]) %>%
                    as.matrix

                ori_intcount_input <- fread(
                    file.path(
                        outDir,
                        folders$int,
                        paste("pivotORI", tissues[i], "Rdelaun.csv", sep = "_")
                )) %>%
                    filter(cell1type != 'dirt') %>%
                    select(colnames(.)[colnames(.) != 'dirt']) %>%
                    as.data.frame %>%
                    column_to_rownames(colnames(.)[1]) %>%
                    as.matrix

                b         <- make_matr(colnames(ori_intcount_input))

                leftcell  <- matrix(
                    apply(
                        expand.grid(colnames(ori_intcount_input),
                                    colnames(ori_intcount_input)),
                        1,
                        function(x){
                            x <- unlist(x)
                            reg_counts[get(annotCol) == x[1],][[tissues[i]]]
                        }) %>% unlist,
                    nrow = ncol(ori_intcount_input)
                )

                rightcell <- matrix(
                    apply(
                        expand.grid(colnames(ori_intcount_input),
                                    colnames(ori_intcount_input)),
                        1,
                        function(x){
                            x <- unlist(x)
                            reg_counts[get(annotCol) == x[2],][[tissues[i]]]
                        }) %>% unlist,
                    nrow = ncol(ori_intcount_input)
                )

                combo1 <- data.table(a = upperTrivec(b),
                                     b = upperTrivec(leftcell),
                                     c = upperTrivec(rightcell),
                                     d = upperTrivec(ori_intcount_input),
                                     e = upperTrivec(lglklyhd_input))

                colnames(combo1) <- c('annotation',
                                      paste(tissues[i],"leftcell_count", sep="_"),
                                      paste(tissues[i],"rightcell_count", sep="_"),
                                      paste(tissues[i],"ori_int_count", sep="_"),
                                      paste(tissues[i],"likelihood_ratio", sep="_"))

                if(nrow(combo) > 0){
                    combo <- combo %>% full_join(combo1, by = 'annotation') %>% as.data.table
                } else{
                    combo <- combo1
                }
            }

            fwrite(combo, file.path(outDir, folders$int, paste0(
                "new_full_pairwise_analysis_UPPERtr_",
                efacts[j],
                ".csv")))
        }
        cat( file = stderr(), 'Ready with new full pairwise analysis ...\n')
        print('Ready with new full pairwise analysis.')
    })

# profiles ------------------------------------------------------------------------------------

    output$profmsg <- renderPrint({

        req(input$data_file)
        outDir      <- path()
        dt          <- dt()
        annotCol    <- ann_colname  #annot_col()
        tissueCol   <- tile_colname  #tissue_col()
        tissues     <- dt[[tissueCol]] %>% unique
        normProfs   <- data.table()
        noClust     <- no_clust()
        ef2reg_map  <- dt[, .(efact, region)] %>% unique %>% .[order(efact)]
        efacts      <- sort(unique(ef2reg_map$efact))
        row_order   <- character()
        for(ef in efacts){
            z <- ef2reg_map[efact == ef, region]
            z <- paste0('reg', sort(as.integer(gsub('reg', '', z))))
            row_order <- c(row_order, z)
        }

        for(fname in dir(file.path(outDir, folders$del), full.names = TRUE)){
            if(grepl('^reg.+_Rdelaun.csv$', basename(fname))){
                #print(fname)
                input <- fread(fname)
                regionid <- strsplit(basename(fname), "_" )[[1]][1]

                input[, cell := apply(input, 1, function(x){
                    x <- unlist(x)
                    paste(regionid, as.integer(x[1]), as.integer(x[2]), sep = '_')
                })]
                input[, dummy := 1]
                profiles <- dcast(input,
                                  cell ~ cell2type,
                                  value.var = 'dummy',
                                  fun.aggregate = length) %>% as.data.table

                profiles_norm <- bind_cols(profiles[, .(cell)], (profiles[, -1]/rowSums(profiles[, -1])))
                fwrite(profiles_norm, file.path(outDir, folders$int, paste0(regionid, "_neighborprofiles.csv")))
                normProfs     <- bind_rows(normProfs, profiles_norm) %>% replace(., is.na(.), 0) %>% as.data.table
            }
        }

        cat(file=stderr(), 'Made cell profiles ...\n')
        print('Made cell profiles.')

        res_km <- kmeans(normProfs[, -'cell'], noClust, iter.max = 100)

        d      <- t(res_km$centers) # celltypes vs clusters
        b      <- data.table(cluster = res_km$cluster, cellID = normProfs$cell)

        dt     <- dt %>% unite('cellID', tissueCol, `X:X`, `Y:Y`, sep = '_', remove = FALSE)

        # add cluster id to ori data
        dt     <- dt %>% left_join(b, by = "cellID") %>% as.data.table %>% .[!is.na(cluster)]
        fwrite(dt, file.path(outDir, 'dt_with_clusters.txt'), sep = '\t')

        quantile.range <- quantile(as.numeric(res_km$centers), probs = seq(0, 1, 0.01))
        #quantile.range <- quantile(as.numeric(log(d + min(d[d>0]))), probs = seq(0, 1, 0.01))
        palette.breaks <- seq(quantile.range["5%"], quantile.range["98%"],
                              (quantile.range["98%"]- quantile.range["5%"])/1000)

        ## use http://colorbrewer2.org/ to find optimal divergent color palette (or set own)
        color.palette  <- colorRampPalette(c("black","blue","green"))(length(palette.breaks) - 1)

        png(file.path(outDir, folders$niche, paste("niches",'png',sep='.')),
            height=1100, width=2000, res = 300, pointsize=10)

        heatmap.2(d, labCol=colnames(d), dendrogram="none", col=color.palette, trace="none",
                  breaks=palette.breaks, margins=c(3,20), cexRow=0.7, cexCol=0.5)
        dev.off()

        # to prevent printing to the plots pane, create a file we don't need
        png(file.path(outDir, folders$niche, paste("waste",'png',sep='.')),
            height=1100, width=2000, res = 300, pointsize=10)
        d        <- d[,order(as.numeric(colnames(d)))]

        nichemap <- heatmap.2(d, labCol=colnames(d), dendrogram="none", col=color.palette,
                              trace="none", breaks=palette.breaks, margins=c(3,20),
                              cexRow=0.7, cexCol=0.5)
        dev.off()
        unlink(file.path(outDir, folders$niche, paste("waste",'png',sep='.')))

        nicheorder <- nichemap$colInd
        cellorder  <- nichemap$rowInd

        d          <- d[,nicheorder]
        d          <- d[cellorder,]
        d          <- d[dim(d)[1]:1,]

        # profiles in niches
        per_region_niches <- dt %>% select(tissueCol, annotCol, nichename = cluster) %>% as.data.table
        per_region_niches[, dummy := 1]

        # cell count per region + celltype for each cluster (columns)
        cells_in_niches   <- dcast(per_region_niches,
                                   as.formula(paste(tissueCol, '+', annotCol, '~', 'nichename')),
                                   value.var = 'dummy',
                                   fun.aggregate = sum)

        setcolorder(cells_in_niches, c(1,2,2+nicheorder))

        cells_in_niches_norm <- bind_cols(cells_in_niches[,1:2],
                                          cells_in_niches[, -c(1:2)] / rowSums(cells_in_niches[, -c(1:2)]))

        fwrite(cells_in_niches, file.path(outDir, folders$int, paste("cells_in_niches.cst",sep = "")))
        fwrite(cells_in_niches_norm, file.path(outDir, folders$int, paste("cells_in_niches_norm.cst",sep = "")))

        celltypes <- unique(cells_in_niches_norm[[annotCol]]) %>% unlist %>% unname


        color.palette <- circlize::colorRamp2(seq(from = 0, to = 1, length.out = 100),
                                              colorRampPalette(c("royalblue4","yellow", 'red'))(100))

        # cell types across clusters by region
        for (index in 1:length(celltypes)){
            #print(index)
            addon <- cells_in_niches_norm[get(annotCol) == celltypes[index], ] %>%
                select(region, colnames(d)) %>%
                column_to_rownames(colnames(.)[1]) %>%
                as.matrix() %>% .[row_order[row_order %in% rownames(.)], ]
            kaugm <- rbind(d, addon) # what is '5'

            df <- data.frame(samples = c(rep(annotCol, dim(d)[1]), ef2reg_map[region %in% rownames(addon), efact]),
                                         stringsAsFactors = FALSE)
            sidebar_col <- c( 'red', 'green')
            if(length(unique(df$samples)) > 2){
                sidebar_col <- brewer.pal(length(unique(df$samples)), 'Set1')
            }

            names(sidebar_col) <- unique(df$samples)

            ha <- HeatmapAnnotation(df = df, which = 'row', col = list(samples = sidebar_col))

            png(file.path(outDir,
                          folders$niche,
                          paste0(clean_path(celltypes[index]),".png")),
                height=6000, width=4500, res = 300, pointsize=10)
            print(
                ha + Heatmap(kaugm/rowSums(kaugm),
                             row_order = rownames(kaugm),
                             column_order = colnames(kaugm),
                             col = color.palette,
                             row_split = factor(df$samples, levels = df$samples %>% unique),
                             cluster_row_slices = F,
                             column_title = paste('Niche composition of regions for :', celltypes[index]),
                             heatmap_legend_param = list(title = 'fraction'))
            )
            dev.off()
        }
        cat(file = stderr(), 'Finished niche clustering ...\n')
        print('Finished niche clustering')
    })


    # ---------------- cell distances --------------------- #############
    output$celldistmsg <- renderPrint({

        req(input$data_file)
        outDir      <- path()
        dt          <- dt()
        annotCol    <- ann_colname  #annot_col()
        celltypes   <- dt[[annotCol]] %>% unique %>% sort
        tissueCol   <- tile_colname  #tissue_col()
        tissues     <- dt[[tissueCol]] %>% unique
        ef2reg_map  <- dt[, .(efact, region)] %>% unique %>% .[order(efact)]
        efacts      <- sort(unique(ef2reg_map$efact))
        convFactor  <- conv_factor()

        min_count   <- 10
        num_bins    <- 5

        dt[, rowid := .I]
        dtd <- data.table()

        # calculate the minimal distance of each cell to cells of every cell type
        for(r in tissues){
            d <- bind_cols(dt[ get(tissueCol) == r, .(rowid)],
                           nndist(dt[get(tissueCol) == r, .(`X:X`, `Y:Y`)],
                                  by = as.factor(dt[get(tissueCol) == r, celltypes])) %>%
                               as.data.table)
            d[[tissueCol]] <- r
            dtd <- bind_rows(dtd, d)
        }

        dt <- dt %>% left_join(dtd, by = c('rowid', 'region')) %>% as.data.table

        # compare dist distributions of celltype1 vs celltype2 in each group
        # 2 different ways
        cs <- data.table()
        sapply(celltypes, function(cl1){
            sapply(celltypes, function(cl2){
                print(paste(cl1, cl2))
                # efacts for this cell type
                efcts <- dt[celltypes == cl1, efact] %>% unique %>% sort
                if(sum(!is.na(dt[celltypes == cl1, get(cl2)])) > min_count){
                    q <- quantile(dt[celltypes == cl1, get(cl2)] , seq(0, 1, (1/num_bins)), na.rm = T)
                    dt[celltypes == cl1, bin := cut(dt[celltypes == cl1, get(cl2)], breaks=q, include.lowest=T)]
                    tab <- with(dt[celltypes == cl1], table(bin, efact))
                    aa <- bind_cols(
                        data.table(
                            cell1 = cl1,
                            cell2 = cl2,
                            chis_pval = suppressWarnings(chisq.test(tab)$p.value)
                        ),
                        dt[celltypes == cl1, median(get(cl2), na.rm = T), by = .(efact, celltypes)] %>%
                            dcast(celltypes ~ efact, value.var = 'V1') %>%
                            .[,-'celltypes'] %>% setnames(., colnames(.), paste('med', colnames(.), sep = '_')),
                        dt[celltypes == cl1, mad(get(cl2), na.rm = T), by = .(efact, celltypes)] %>%
                            dcast(celltypes ~ efact, value.var = 'V1') %>%
                            .[,-'celltypes'] %>% setnames(., colnames(.), paste('mad', colnames(.), sep = '_')),
                        dt[celltypes == cl1 & !is.na(get(cl2)), .N, by = .(efact, celltypes)] %>%
                            dcast(celltypes ~ efact, value.var = 'N') %>%
                            .[,-'celltypes'] %>% setnames(., colnames(.), paste('count', colnames(.), sep = '_'))
                    )
                    dt[, bin := NULL]

                    # print('past aa...')
                    tt <- data.table()
                    if(length(efcts) >= 2){
                        sapply(1:(length(efcts)-1), function(g1){
                            sapply((g1+1):length(efcts), function(g2){
                                # print(c(g1, '   ', g2))
                                g1 <- efcts[g1]
                                g2 <- efcts[g2]
                                # print(paste(g1, g2))
                                tt <<- bind_rows(tt,
                                                 data.table(
                                                     cell1 = cl1,
                                                     cell2 = cl2,
                                                     comp  = paste(g1, 'vs', g2, 'mw', sep = '_'),
                                                     pval  = ifelse(
                                                         sum(!is.na(dt[celltypes == cl1 & efact == g1, get(cl2)])) > min_count &&
                                                             sum(!is.na(dt[celltypes == cl1 & efact == g2, get(cl2)])) > min_count,
                                                         suppressWarnings(
                                                             wilcox.test(dt[celltypes == cl1 & efact == g1, get(cl2)],
                                                                         dt[celltypes == cl1 & efact == g2, get(cl2)])$p.value
                                                         ),
                                                         NA_real_
                                                     )
                                                 ))
                            })
                        })
                        tt <- tt %>% dcast(cell1 + cell2 ~ comp, value.var = 'pval')
                        aa <- aa %>% left_join(tt, by = c('cell1', 'cell2')) %>% as.data.table()

                        # print('past tt')
                        dd <- data.table()
                        sapply(1:(length(efcts)-1), function(g1){
                            sapply((g1+1):length(efcts), function(g2){
                                g1 <- efcts[g1]
                                g2 <- efcts[g2]
                                # print(paste(g1, g2))
                                dd <<- bind_rows(dd,
                                                 data.table(
                                                     cell1 = cl1,
                                                     cell2 = cl2,
                                                     comp  = paste(g1, 'vs', g2, 'l2r', sep = '_'),
                                                     l2r   = ifelse(
                                                         !is.na(aa[, paste0('med_', g1), with = F])
                                                         && !is.na(aa[, paste0('med_', g2), with = F]),
                                                         log2(aa[, paste0('med_', g1), with = F] / aa[, paste0('med_', g2), with = F]),
                                                         NA_real_
                                                     ) %>% unlist %>% unname
                                                 ))
                                #print(dd)
                            })
                        })

                        dd <- dd %>% dcast(cell1 + cell2 ~ comp, value.var = 'l2r')
                        aa <- aa %>% left_join(dd, by = c('cell1', 'cell2')) %>% as.data.table()
                        # print('past dd')
                    }

                } else {
                    aa <- data.table(
                        cell1 = cl1, cell2 = cl2, chis_pval = NA
                    )
                }
                cs <<- bind_rows(cs, aa)
            })
        })

        # fdr correction:
        cs <- bind_cols(cs[, -which(grepl('_mw$', colnames(cs))), with = F],
                        cs[, which(grepl('_mw$', colnames(cs))), with = F] %>%
                            unlist %>%
                            p.adjust(method = 'fdr') %>%
                            matrix(ncol = length(which(grepl('_mw$', colnames(cs))))) %>%
                            as.data.table %>%
                            setnames(colnames(cs)[which(grepl('_mw$', colnames(cs)))])
        )

        ks <- data.table()
        sapply(celltypes, function(cl1){
            sapply(celltypes, function(cl2){
                sapply(1:(length(efacts)-1), function(n1){
                    sapply((n1+1):length(efacts), function(n2){
                        g1 <- efacts[n1]
                        g2 <- efacts[n2]
                        print(paste(cl1, cl2, n1, n2))
                        print(paste(cl1, cl2, g1, g2))
                        if(sum(!is.na(dt[celltypes == cl1 & efact == g1, get(cl2)])) > 10 &
                           sum(!is.na(dt[celltypes == cl1 & efact == g2, get(cl2)])) > 10){
                            d1 <- dt[celltypes == cl1 & efact == g1, get(cl2)]
                            d2 <- dt[celltypes == cl1 & efact == g2, get(cl2)]
                            ks <<- bind_rows(ks, data.table(
                                cell1 = cl1, cell2 = cl2, group1 = g1, group2 = g2,
                                pval = suppressWarnings(ks.test(d1, d2)$p.value)
                            ))
                        } else {
                            ks <<- bind_rows(ks, data.table(
                                cell1 = cl1, cell2 = cl2, group1 = g1, group2 = g2, pval = NA
                            ))
                        }
                    })
                })
            })
        })

        # ks[, comp := apply(.SD, 1, function(x){paste0(x[1], '_vs_', x[2], '_ks')}), .SDcols = c('group1', 'group2')]
        #ks[, comp := apply(ks[, c('group1', 'group2'), with = F], 1, function(x){paste0(x[1], '_vs_', x[2], '_ks')})]
        ks[, comp := paste0(group1, '_vs_', group2, '_ks')]
        ks[, padj := p.adjust(pval)]
        ks <- ks %>% dcast(cell1 + cell2 ~ comp, value.var = 'padj')

        # result of the comparison
        res <- cs %>%
            left_join(ks,  by = c('cell1', 'cell2')) %>%
            mutate(chis_padj = p.adjust(chis_pval)) %>%
            as.data.table

        resplot <- function(zz, cl1, cl2, ptable = NA, pval = NA, logtr = TRUE, dens = FALSE){

            dpoints <- zz[celltypes == cl1, cl2, with = F] %>% unlist %>% unname

            suppressWarnings(
                if(logtr){
                    p <- ggplot(zz[celltypes == cl1, ], aes(log10(get(cl2)), color = as.factor(efact)))
                    dpoints <- sum(!is.infinite(dpoints) & !is.na(dpoints))
                } else {
                    p <- ggplot(zz[celltypes == cl1, ], aes(get(cl2), color = as.factor(efact)))
                    dpoints <- sum(!is.infinite(dpoints) & !is.na(dpoints))
                }
            )
            tt <- ifelse(dens, 'density', 'histogram')
            title <- paste0(cl1, ' vs. ', cl2, ', ', tt, ' (no of datapoints: ', dpoints, ')')
            if( !is.na(pval)){
                title <- paste0(title, '; pval = ', formatC(pval, format = 'e', digits = 2))
            }
            if(dens){
                p <- p + geom_density(adjust = 3/4)
            } else{
                p <- p  + geom_histogram(bins = 100, position = 'dodge', alpha = 0.5)
            }
            p <- p + xlab(paste0('log10(dist to ', cl2, ')')) +
                labs(color = 'efact') +
                ggtitle(title) +
                theme(plot.title = element_text(hjust = 0.5))
            if('data.frame' %in% class(ptable)){
                p <- p + annotate(geom = "table", x = 0, y = 1, label = list(ptable))
            }
            p
        }

        # draw histograms or density plots of the distributions
        sapply(celltypes, function(c1){
            sapply(celltypes, function(c2){
                pvals_table <- res[cell1 == c1 & cell2 == c2, which(grepl('_mw$', colnames(res))), with = F] %>%
                    as.data.frame() %>%
                    t %>%
                    as.data.table(keep.rownames = T) %>%
                    .[, .(comparison = rn, pval = V1)]

                gr <- resplot(dt, c1, c2, ptable = pvals_table, dens = dens_his() == 'dens')
                pdf(file = file.path(outDir,
                                     folders$dist,
                                     clean_path(paste0('celldist_', c1, '_vs_', c2, '.pdf'))))
                print(gr)
                dev.off()
            })
        })

        fwrite(res, file.path(outDir, folders$dist, 'celldist.txt'), sep = '\t')

        dt <- dt[, c(tissueCol, x_coord_colname, y_coord_colname, 'efact', annotCol, celltypes), with = F]
        colnames(dt)[2:5] <- c(x_col(), y_col(), efact_col(), annot_col())

        # restore previously removed rows
        dt_ori <- fread(inputf_ori()) %>%
            .[, c(annot_col(), tissue_col(), efact_col(), pat_col(), x_col(), y_col()), with = F]

        if(!all(grepl('^reg', dt_ori[[tissue_col()]]))){
            dt_ori[, (tissueCol) := paste0('reg', get(tissue_col()))]
        } else {
            dt_ori[, (tissueCol) := get(tissue_col())]
        }

        dt_ori[, (tissueCol) := lapply(.SD, function(x){gsub('reg0*([1-9]\\d*)$', 'reg\\1', x)}), .SDcols = tissueCol]

        # celltypes need to be updated in case we filtered out some of them
        # earlier, because of min cell count
        celltypes <- sort(unique(dt_ori[[annot_col()]]))

        dt_ori <- dt_ori %>%
            left_join(dt, by = c(tissueCol, x_col(), y_col(), efact_col(), annot_col())) %>%
            as.data.table %>%
            .[, -c(tissueCol), with = F]

        # if some celltypes were filtered out, we need to add back an emtpy column for that type
        missing <- celltypes[which(!celltypes %in% colnames(dt_ori))]
        if(length(missing) > 0){
            for(k in missing){
                dt_ori[[k]] <- NA
            }
        }
        dt_ori <- setcolorder(dt_ori, c(tissue_col(), x_col(), y_col(), efact_col(),
                                        annot_col(), pat_col(), celltypes))

        fwrite(dt_ori, file.path(outDir, 'dt_celldist_pixel.txt'), sep = '\t')

        # convert distances to microm
        dt_ori_um <- bind_cols(dt_ori[,1:6], dt_ori[, lapply(.SD, function(x){x*convFactor}), .SDcols = 7:dim(dt_ori)[2]])
        fwrite(dt_ori_um, file.path(outDir, 'dt_celldist_um.txt'), sep = '\t')

        set.seed(12345)

        # earlier version of dist_sign...
        dist_sign1 <- function(dtd, ct2, n_mean = 1000, ps = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)){

            # reshape data table and keep only rows where cell2type == ct2
            dtdf <- dtd %>%
                melt(id.vars = c(annot_col(), tissue_col(), efact_col(), x_col(), y_col(), pat_col()),
                     variable.name = 'cell2type',
                     value.name = 'distance') %>%
                .[!is.na(distance)] %>%
                .[cell2type == ct2]

            res <- data.table()
            sapply(celltypes, function(ct1){
                cat(file = stderr(), paste(ct1, ct2, sep = ' --- '))
                res <<- bind_rows(res, dtdf[, # sample distances to non-ct1 cells as many times as the
                              # number of ct1 cells, get their mean, n_mean times,
                              # this is all by region
                            .(dist = replicate(n_mean,
                                               mean(
                                                   sample(.SD[get(annot_col()) != ct1, distance],
                                                          nrow(.SD[get(annot_col()) == ct1]),
                                                          replace = T)),
                                               simplify = T),
                              N_all = .N), # number of distances to ct2 per region
                            c(tissue_col())][order(get(tissue_col()))] %>%
                    full_join(dtdf[get(annot_col()) == ct1,
                                   # ave distance between ct1 and ct2 and number of distances
                                   # between them per region
                                   .(mean_ct1 = mean(distance),
                                     N_ct1 = .N), c(tissue_col())],
                              by = tissue_col()) %>%
                    as.data.table %>%
                    .[, # calculate various statistics
                      .(N_ct1 = unique(N_ct1),
                        mean_ct1 = unique(mean_ct1),
                        quant = sum(mean_ct1 > dist, na.rm = T),
                        d_mean = mean(dist, na.rm = T), # mean of n_mean random means
                        d_sd = sd(dist, na.rm = T), # sd of n_mean random means
                        quants = quantile(dist, probs = ps, na.rm = T), # the quantile of n_mean random means
                        cols = paste0('d_', ps)), # which quantile
                      c(tissue_col(), 'N_all')] %>%
                    # reformat table
                    dcast(formula(paste(tissue_col(), '+ mean_ct1 + N_ct1 + N_all + quant + d_mean + d_sd ~ cols')),
                          value.var = 'quants') %>%
                    # calculate z-score of actual distance to random distribution
                    .[, zscore := (mean_ct1 - d_mean) / d_sd] %>%

                    full_join( # add back columns that we need in future
                        dtdf[, -c('distance', x_col(), y_col(), annot_col()), with = F] %>%
                            unique() %>%
                            mutate(!!annot_col() := ct1),
                        by = tissue_col()) %>%
                    as.data.table %>%
                    setcolorder(c(annot_col(), 'cell2type', tissue_col(), pat_col(), efact_col(),
                                  'N_ct1', 'N_all', 'quant', 'zscore', 'mean_ct1')))
            })
            invisible(res)
        }

        ac <- annot_col()
        tc <- tissue_col()
        xc <- x_col()
        yc <- y_col()
        ec <- efact_col()
        pc <- pat_col()

        dist_sign <- function(dtd, ct2, by_what, n_mean = 1000, ps = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)){
            # calculate significance of ave min distance between ct1s and ct2 per region
            # based on the ave distances of randonm non-ct1 distances to ct2

            # reshape data table and keep only rows where cell2type == ct2

            cols_to_keep <- c(ac, tc, ec, xc, yc, pc)
            col_order <- c(ac, 'cell2type', tc, pc, ec, 'N_ct1', 'N_all', 'quant', 'zscore', 'mean_ct1')
            if(by_what == ec){
                cols_to_keep <- c(ac, ec, xc, yc)
                col_order <- c(ac, 'cell2type', ec, 'N_ct1', 'N_all', 'quant', 'zscore', 'mean_ct1')
            }

            dtdf <- dtd %>%
                melt(id.vars = cols_to_keep,
                     variable.name = 'cell2type',
                     value.name = 'distance') %>%
                .[!is.na(distance)] %>%
                .[cell2type == ct2]

            res <- data.table()
            if(nrow(dtdf) > 0){
                #sapply(celltypes, function(ct1){
                sapply(sort(unique(dtdf[[ac]])), function(ct1){
                    cat(file = stderr(), paste(ct1, ct2, sep = ' --- '), '\n')
                    res <<- bind_rows(res, dtdf[, # sample distances to non-ct1 cells as many times as the
                                                # number of ct1 cells, get their mean, n_mean times,
                                                # this is all by region
                                                .(dist = replicate(n_mean,
                                                                   mean(
                                                                       sample(.SD[get(ac) != ct1, distance],
                                                                              nrow(.SD[get(ac) == ct1]),
                                                                              replace = T)),
                                                                   simplify = T),
                                                  N_all = .N), # number of distances to ct2 per region
                                                c(by_what)][order(get(by_what))] %>%
                                          full_join(dtdf[get(ac) == ct1,
                                                         # ave distance between ct1 and ct2 and number of distances
                                                         # between them per region
                                                         .(mean_ct1 = mean(distance),
                                                           N_ct1 = .N), c(by_what)],
                                                    by = by_what) %>%
                                          as.data.table %>%
                                          .[, # calculate various statistics
                                            .(N_ct1 = unique(N_ct1),
                                              mean_ct1 = unique(mean_ct1),
                                              quant = sum(mean_ct1 > dist, na.rm = T),
                                              d_mean = mean(dist, na.rm = T), # mean of n_mean random means
                                              d_sd = sd(dist, na.rm = T), # sd of n_mean random means
                                              quants = quantile(dist, probs = ps, na.rm = T), # the quantile of n_mean random means
                                              cols = paste0('d_', ps)), # which quantile
                                            c(by_what, 'N_all')] %>%
                                          # reformat table
                                          dcast(formula(paste(by_what, '+ mean_ct1 + N_ct1 + N_all + quant + d_mean + d_sd ~ cols')),
                                                value.var = 'quants') %>%
                                          # calculate z-score of actual distance to random distribution
                                          .[, zscore := (mean_ct1 - d_mean) / d_sd] %>%

                                          full_join( # add back columns that we need in future
                                              dtdf[, -c('distance', xc, yc, ac), with = F] %>%
                                                  unique() %>%
                                                  mutate(!!ac := ct1),
                                              by = by_what) %>%
                                          as.data.table %>%
                                          setcolorder(col_order))
                })
            }
            invisible(res)
        }


        cores <- detectCores()
        cl    <- makeCluster(min(c(length(celltypes), cores)),
                             outfile = file.path(outDir, folders$int, paste0('dopar1_log')))
        registerDoParallel(cl)
        set.seed(12345)
        dist_signif <- foreach(
            i = 1:length(celltypes),
            .combine = 'bind_rows',
            # .export = c('dist_sign'),
            .packages = c('data.table', 'dplyr')
        ) %dopar% {
            dist_sign(dt_ori_um, celltypes[i], tc)
        }
        fwrite(dist_signif, file.path(outDir, paste0('dt_celldist_um_signif_by', tc, '.txt')), sep = '\t')

        set.seed(12345)
        dist_signif <- foreach(
            i = 1:length(celltypes),
            .combine = 'bind_rows',
            # .export = c('dist_sign'),
            .packages = c('data.table', 'dplyr')
        ) %dopar% {
            dist_sign(dt_ori_um, celltypes[i], ec)
        }
        fwrite(dist_signif, file.path(outDir, paste0('dt_celldist_um_signif_by', ec, '.txt')), sep = '\t')
        stopCluster(cl)

        random_ratios1 <- function(dt, col, tue, n_repl = 100, ps = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)){

            # significance of right_cell - middle_cell distance / middle_cell - left_cell distance, per region
            # The random samples are drawn from all right_cell / left_cell ratios in the region

            ulabels    <- unique(sort(dt[[ac]]))
            fullres_dt <- data.table()

            cols_include <- c(tc, ec, pc)
            if(col == ec){
                cols_include <- col
            }

            for(left_cell in ulabels){
                for(right_cell in ulabels){
                    print(paste('Method1', col, tue, left_cell,  right_cell, sep = ' -- '))
                    # skip when left == right, ratio would be 1
                    if(left_cell == right_cell) next

                    dt1 <- dt[, c(ac, left_cell, right_cell), with = F]
                    dt1 <- dt1[, ratio := dt1[[right_cell]] / dt1[[left_cell]]] %>%
                        .[, c(ac, 'ratio'), with = F]
                    #%>%
                    #    .[!is.na(ratio)]
                    for(middle_cell in ulabels){
                        print(paste(col, tue, left_cell, middle_cell, right_cell, sep = ' -- '))
                        sample_size <- nrow(dt[get(ac) == middle_cell,])

                        randrats <- replicate(n_repl, sample(dt1[, ratio], size = sample_size) %>% mean(na.rm = T))

                        mean_val <- dt1[get(ac) == middle_cell, ratio] %>% mean(na.rm = T)

                        res_dt <- bind_cols(unique(dt[, cols_include, with = F]),
                                            data.table(cell_left = left_cell, cell_mid = middle_cell,
                                                       cell_right = right_cell, N_mid = sample_size,
                                                       N_tot = nrow(dt1), ratio_mean = mean_val,
                                                       quant = sum(mean_val > randrats, na.rm = T),
                                                       mean_rand = mean(randrats),
                                                       sd_rand = sd(randrats)),
                                            setnames(
                                                data.table(t(quantile(randrats, ps, na.rm = T))),
                                                paste('d', ps, sep = '_')
                                            ))

                        fullres_dt <- bind_rows(fullres_dt, res_dt)
                    }
                }
            }
            invisible(fullres_dt)
        }

        random_ratios2 <- function(dt, col, tue, n_repl = 100, ps = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)){

            # another way of determining the significance of
            # right_cell - middle_cell distance / middle_cell - left_cell distance, per region
            # the random distances are drawn from the regions distance matrix, where the
            # row and column labels are shuffled

            # dist_mat   <- dist(dt[, c(xc, yc), with = F]) %>% as.matrix() # dist of each cell to each cell in region
            labels     <- dt[[ac]]
            ulabels    <- unique(sort(labels))
            fullres_dt <- data.table()

            # don't want 0s in this matrix
            # diag(dist_mat) <- max(dist_mat)

            mean_val1 <- function(dm, labs, lc, rc, mc){
                # nonzero_min <- function(x){ min(x[x>0])}
                mix <- which(labs == mc)
                rix <- which(labs == rc)
                lix <- which(labs == lc)
                mean(rowMins(dm[mix, rix], value = T) /
                         rowMins(dm[mix, lix], value = T), na.rm = T)

                #(apply(dm[mix, rix, drop = F], 1, nonzero_min) /
                #     apply(dm[mix, lix, drop = F], 1, nonzero_min)) %>%
                #     mean(na.rm = T)
            }


            mean_val <- function(dm, labs, lc, rc, mc){
                mix <- which(labs == mc)
                rix <- which(labs == rc)
                lix <- which(labs == lc)
                mean(
                    rowMins(dista(dm[mix], dm[rix]), value = T)
                    /
                        rowMins(dista(dm[mix], dm[lix]), value = T),
                    na.rm = T
                )
            }

            for(left_cell in ulabels){
                for(right_cell in ulabels){

                    # skip when left == right, ratio would be 1
                    if(left_cell == right_cell) next

                    for(middle_cell in ulabels){

                        ctypes             <- c(left_cell, middle_cell, right_cell)
                        sample_size        <- sum(labels == middle_cell)
                        mean_min_dist      <- NA
                        mean_rand_min_dist <- NA
                        print(paste(c('Method2', col, tue, ctypes, memory.size()), collapse = ' -- '))

                        if(all(ctypes %in% labels)){
                            # mean_min_dist <- mean_val(dist_mat, labels, left_cell, right_cell, middle_cell)
                            mean_min_dist <- mean_val(dt[, c(xc, yc), with = F], labels, left_cell, right_cell, middle_cell)
                            # mean_rand_min_dist <- replicate(n_repl,
                            #                                 mean_val(dist_mat, sample(labels),
                            #                                          left_cell, right_cell, middle_cell),
                            #                                 simplify = 'vector')
                            mean_rand_min_dist <- replicate(n_repl,
                                                            mean_val(dt[, c(xc, yc), with = F], sample(labels),
                                                                     left_cell, right_cell, middle_cell),
                                                            simplify = 'vector')

                            res_dt <- bind_cols(unique(dt[, c(tc, ec, pc), with = F]),
                                                data.table(cell_left = left_cell, cell_mid = middle_cell,
                                                           cell_right = right_cell, N_mid = sample_size,
                                                           N_tot = length(labels), ratio_mean = mean_min_dist,
                                                           quant = sum(mean_min_dist > mean_rand_min_dist, na.rm = T),
                                                           mean_rand = mean(mean_rand_min_dist),
                                                           sd_rand = sd(mean_rand_min_dist)),
                                                setnames(
                                                    data.table(t(quantile(mean_rand_min_dist, ps, na.rm = T))),
                                                    paste('d', ps, sep = '_')
                                                ))

                            fullres_dt <- bind_rows(fullres_dt, res_dt)

                        } else {
                            cat(file = stderr(), paste(c(tue, 'missing' ,ctypes[which(!ctypes %in% labels)]), collapse = ' -- '), '\n')
                        }
                    }
                }
            }
            invisible(fullres_dt)
        }

        spots  <- unique(dt_ori_um[[tc]]) %>% sort
        efacts <- unique(dt_ori_um[[ec]]) %>% sort
        cl     <- makeCluster(min(c(length(spots), cores)),
                              outfile = file.path(outDir, folders$int, 'dopar2_log'))
        registerDoParallel(cl)

        cat(file = stderr(), 'Calculating ratios by spot, method1.\n')
        set.seed(12345)
        rands <- foreach(
            i = 1:length(spots),
            .combine = 'bind_rows',
            .packages = c('data.table', 'dplyr')
        ) %dopar% {
            random_ratios1(dt_ori_um[get(tc) == spots[i], ], tc, spots[i])
        }
        fwrite(rands, file.path(outDir, 'dt_celldist_ratios_bySpot_meth1.txt'), sep = '\t')

        cat(file = stderr(), 'Calculating ratios by groups, method1.\n')
        set.seed(12345)
        rands <- foreach(
            i = 1:length(efacts),
            .combine = 'bind_rows',
            .packages = c('data.table', 'dplyr')
        ) %dopar% {
            random_ratios1(dt_ori_um[get(ec) == efacts[i], ], ec, efacts[i])
        }
        fwrite(rands, file.path(outDir, 'dt_celldist_ratios_byGroup_meth1.txt'), sep = '\t')

        stopCluster(cl)

        cl     <- makeCluster(min(c(length(spots), cores/2)),
                              outfile = file.path(outDir, folders$int, 'dopar3_log'))
        registerDoParallel(cl)

        cat(file = stderr(), 'Calculating ratios by spot, method2.\n')
        set.seed(12345)
        rands <- foreach(
            i = 1:length(spots),
            .combine = 'bind_rows',
            .packages = c('data.table', 'dplyr', 'Rfast')
        ) %dopar% {
            random_ratios2(dt_ori_um[get(tc) == spots[i], ], tc, spots[i])
        }
        fwrite(rands, file.path(outDir, 'dt_celldist_ratios_meth2.txt'), sep = '\t')

        stopCluster(cl)

        cat(file = stderr(), 'Finished cell distance calculations.\n')
        print('Finished cell distance calculations')
    })
}

# Run the application
shinyApp(ui = ui, server = server)
