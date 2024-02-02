library(shiny)
library(shinyhelper) 
library(shinythemes)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(tableHTML)
library(cowplot)


ui <- div(
  # https://stackoverflow.com/questions/70655418/shiny-how-to-make-the-title-of-navbarpage-fill-the-window
  # https://community.rstudio.com/t/is-there-a-way-to-make-the-title-in-navbarpage-function-becomes-a-tab/70876/2
  tags$style(HTML("
        body > div > .container-fluid:nth-of-type(1) {
            margin: 0 auto;
            padding-top: 70px;
        }
        body > div .container-fluid {
            padding-left:10%;
            padding-right:15%;
        }
        body > div > nav .nav.navbar-nav {
            float: right;
            padding-right:0%;
            padding-left:10%;
            height: 25px;
            min-height:25px !important;
            font-size: 20px !important;
        }
        body > div > nav .navbar-brand { font-size: 32px; height: 25px; }
        body > div > nav .navbar-nav li a {font-size: 20px; font-weight: bold; }
        .footer {position:fixed; left:0; bottom:0; padding-left:10%; width:100%; background-color:#f8f8f8;
                       height:5%; padding-top:8px; border-top:1px solid #e7e7e7; text-align:left;}
        ")),
  
  navbarPage(
    position = "fixed-top",
    fluid = TRUE,
    id = "homepage",
    title = actionLink("title", "ImParalog"),
    selected = "ParalogSSEA",
    
    footer = tags$div(class ="footer", p(HTML(paste0("Please cite: ", a("Dong, C., Zhang, F., He, E. et al. In-silico Prediction of Synergistic Paralogs Enhancing Cancer Immunotherapy. JournalName (2024)", href = "wwww.google.com", target="_blank"), ".")), style="color:black;font-size:18px;word-spacing:0;")),
    
    tabPanel("ParalogSSEA",
             tabsetPanel(
               tabPanel("Summary",
                        h4("The paralog sgRNA set enrichment analysis (pSSEA) methods thoughts"),
                        br(),
                        fluidRow(column(1),
                                 column(11,
                                        div(img(src='pic_1.png', height="400px", width="auto"), style="text-align: left;"),
                                 )),
                        br(),
                        br(),
                        ),
               tabPanel(title = "Query this paper", 
                        br(),
                        fluidRow(
                          column(3,
                                 fluidRow(
                                   column(6, style=list("padding-right: 5px;"),
                                          selectInput(inputId = "pssea_gene1", label = "Gene 1:", "")),
                                   column(6, style=list("padding-left: 10px;"),
                                          selectInput(inputId = "pssea_gene2", label = "Gene 2:", ""))),
                                 fluidRow(
                                   column(6, style=list("padding-right: 5px;"),
                                          selectInput(inputId = "pssea_cell", label = "Cell model:", "")),
                                   column(6)),
                                 br()),
                          column(9,
                                 hr(),
                                 fluidRow(
                                   column(6, 
                                          uiOutput("pssea_heatmap_plot.ui"),
                                   ),
                                   column(6, textOutput("pssea_deseq_txt"), #verbatimTextOutput
                                          #tableHTML_output("pssea_deseq_tab"
                                          DT::dataTableOutput('pssea_deseq_tab', width = "auto", height = "auto")
                                   )
                                 ),
                                 hr(),
                                 fluidRow(
                                   column(6, 
                                          uiOutput("pssea_rank_plot.ui"),
                                   ),
                                   column(6,
                                          DT::dataTableOutput('pssea_tab', width = "auto", height = "auto"),
                                   )
                                 ),
                          ))
               ),
               tabPanel(title = "Run with your data", "This section is coming soon.")
             )),
    
    tabPanel(title = "Paralog Prediction",
             br(),
             fluidRow(
               column(3,
                      fluidRow(
                        column(6, style=list("padding-right: 5px;"),
                               selectInput(inputId = "pred_gene1", label = "Gene 1:", "")),
                        column(6, style=list("padding-left: 10px;"),
                               selectInput(inputId = "pred_gene2", label = "Gene 2:", ""))),
                      br(),
                      h5("Download ALL predictions", style="color:black"),
                      downloadButton("pred_downloadData", label = "Download", class="download_this")),
               column(9,
                      fluidRow(column(12, align="center", 
                                      hr(),
                                      uiOutput("pred_rank_plot.ui"),
                                      hr(),
                                      fluidRow(
                                        column(7, 
                                               uiOutput("pred_feat_shap_plot.ui")
                                        ),
                                        column(5,
                                               DT::dataTableOutput('pred_tab', width = "80%", height = "auto"),
                                        )),
                                 ))
               )),
    ),
    
    tabPanel(
      title = "Project Description",
      br(),
      p(HTML("The emerging immunotherapy, especially immune checkpoint blockade (ICB) and chimeric antigen receptor T cell therapy (CAR-T), has revolutionized cancer treatments with remarkable responses in advanced tumors. Still, a large fraction of patients fails to respond to immunotherapy. Recent effort has been made in using CRISPR knockout or activating to screening target to boost T cell effector function and further leverage the immune killing function. However, manipulating a single gene might still be hard to overcome the resistance due to the compromise effect of other genes. Paralogs, derived from the same ancestors and enjoy similar functionality, are reported with synthetic lethal interactions. Simultaneously targeting these paralogs can act synergistically to activate the T cell's function. However, a systematic analysis of paralog pairs that can enhance cancer immunotherapy has yet to be conducted. In this study, we developed a novel machine learning model to predict potential paralog pairs capable of enhancing cancer immunotherapy by incorporating features from aspects as gene characteristics, sequence and structure similarity, protein–protein interaction (PPI) networks, and gene co-evolution information. Firstly, true paralog pairs were identified from immunotherapy pool screen studies using an adopted Kolmogorov-Smirnov test (KS) method. We demonstrated our method enables the identification of paired paralogs beyond traditional single-gene ranking methods and validated through experimental approaches. Subsequently, we constructed an ensemble learning XGBoost classifier to predict true positive immunotherapy paralog pairs. To further validate our findings, we tested the top predicted paralog pairs using external datasets and experimental approaches. The outcomes of this study are expected to provide valuable insights for prioritizing novel targets and inspire effective combination therapeutic strategies for cancer treatment."), style="color:black;font-size:18px;"),
      br(),
      #img(src='pic_2.png', height="50%", width="55%", text-align="center"), #text-align: center;
      fluidRow(column(2),
               column(10,
                      div(img(src='pic_2.png', height="400px", width="auto"), style="text-align: left;"),
               )),
      br(),
      p(HTML(paste0("Please cite: ", a("GEO", href = "wwww.google.com", target="_blank"), ".")),
        style="color:black;font-size:18px;")),
    
    tabPanel(
      title = "About",
      h4("For additional information, please contact Sidi Chen (sidi.chen[at]yale.edu)"),
      br(),
      h4("ShinyApp developed and maintained by Chuanpeng Dong (chuanpeng.dong[at]yale.edu)")
      )
  )
)


pred = readRDS("XGB_predict.rds")
genelist = readRDS("genelist.rds")
tko_genedata<- readRDS('TKO_queryable_genedata.rds')
tko_deseq<- readRDS('TKO_DESeq_result.rds')
tko_pssea<- readRDS('TKO_pSSEA_result.rds')
tko_rawdata<- readRDS('TKO_rawdata_logRPM.rds')


server <- function(input, output, session){
  
  observeEvent(input$title, {
    updateNavbarPage(session, "homepage", "Project Description")
  })
  
  observe_helpers() 
  optCrt="{ option_create: function(data,escape) {return('<div class=\"create\"><strong>' + '</strong></div>');} }" 
  updateSelectizeInput(session, inputId="pssea_gene1", choices = genelist$Symbol, server = TRUE, 
                       selected = 'RAPGEF1', options = list(maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, inputId="pssea_gene2", choices = genelist$Symbol, server = TRUE, 
                       selected = 'RAPGEF2', options = list(maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, inputId="pssea_cell", choices = c('Melanoma - B16',
                                                                  'Breast - 4T1', 'Breast - EMT6',
                                                                  'Colon - CT26', 'Colon - MC38',
                                                                  'Renal - Renca'), server = TRUE, 
                       selected = 'Melanoma - B16', options = list(create = TRUE, persist = TRUE, render = I(optCrt)))
  
  updateSelectizeInput(session, inputId="pred_gene1", choices = genelist$Symbol, server = TRUE, 
                       selected = 'RAPGEF1', options = list(maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, inputId="pred_gene2", choices = genelist$Symbol, server = TRUE, 
                       selected = 'RAPGEF2', options = list(maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
    
  ##### pSSEA sesseion
  pssea_text <- reactiveVal('Please check your input!')
  output$pssea_txt1 <- renderText({ pssea_text() })
  
  pssea_outs<- reactive({ pssea_visz_func(input$pssea_gene1, input$pssea_gene2, input$pssea_cell) })
  
  #output$pssea_deseq_tab<- render_tableHTML({ tableHTML(pssea_outs()[[3]], rownames=FALSE) })
  output$pssea_deseq_tab<- DT::renderDataTable({ DT::datatable( pssea_outs()[[3]], escape=FALSE,rownames= FALSE,
                               options = list(searching=FALSE, pageLength=4,
                                              lengthChange = FALSE)) }) 
  
  pssea_deseq_text<- reactive({
    if (length(pssea_outs()) >2 ) {
      return("DESeq2 raw output")
    } else {
      return(NULL)
    }
  })
  #pssea_deseq_text<- reactiveVal({ ifelse(length(pssea_outs())>2, "DESeq2 raw output", "") })
  output$pssea_deseq_txt <- renderText({ pssea_deseq_text() })
  
  output$pssea_heatmap_plot <- renderPlot({ pssea_outs()[[4]] }) 
  output$pssea_heatmap_plot.ui <- renderUI({ plotOutput("pssea_heatmap_plot",height=200, width="auto") })
  
  output$pssea_tab<- DT::renderDataTable({ DT::datatable( pssea_outs()[[6]], escape=FALSE,rownames= FALSE,
                                                          options = list(searching=FALSE, lengthChange = FALSE)) }) 

  output$pssea_rank_plot <- renderPlot({ pssea_outs()[[5]] }) 
  output$pssea_rank_plot.ui <- renderUI({ plotOutput("pssea_rank_plot",height=360, width="auto") })
  
  
  ##### Prediction sesseion
  pred_text <- reactiveVal('Please check your input!')
  output$pred_txt1 <- renderText({ pred_text() })
  
  pred_outs<- reactive({ predict_visz_func(input$pred_gene1, input$pred_gene2) })
  
  output$pred_rank_plot <- renderPlot({ pred_outs()[[1]] }) 
  output$pred_rank_plot.ui <- renderUI({ plotOutput("pred_rank_plot",height=300, width = "70%") })
  
  output$pred_tab<- DT::renderDataTable({ DT::datatable( pred_outs()[[3]], escape=FALSE,rownames= FALSE,
                                            options = list(searching=FALSE,lengthChange = FALSE)) })
  
  output$pred_feat_shap_plot <- renderPlot({ pred_outs()[[2]] }) 
  output$pred_feat_shap_plot.ui <- renderUI({ plotOutput("pred_feat_shap_plot",height=600, width = "auto" ) })
  
  # Downloadable csv of All predictions
  output$pred_downloadData <- downloadHandler(
    filename = function() { paste("XGBoost_prediction_result.csv", sep ="_") },
    content = function(file) {
      contentt<- pred[,c(2:39)]
      write.csv(contentt, file, row.names = FALSE)
    })
  ############## Func Predict Visz End ###############
  
}


# Function for visualize PSSEA result query
pssea_visz_func<- function(gene1, gene2, cell_sele){
  #print(gene1)
  #print(gene2)
  cells<-  c('B16', '4T1','EMT6','CT26','MC38','Renca');
  names(cells)<- c('Melanoma - B16', 'Breast - 4T1', 'Breast - EMT6',  'Colon - CT26', 'Colon - MC38', 'Renal - Renca')
  cell = cells[cell_sele]
  #cell='B16'
  #gene1='KLHL7'
  #gene2='KLHL13'
  g1_info<- tko_genedata[which(tko_genedata$hgnc_symbol == gene1),]
  g2_info<- tko_genedata[which(tko_genedata$hgnc_symbol == gene2),]
  
  tko_deseq2<- tko_deseq[which(tko_deseq$cell==cell),]
  tko_pssea2<- tko_pssea[which(tko_pssea$cell==cell),]
  id<- paste(min(gene1,gene2), max(gene1,gene2), sep='/')
  #print(id)
  
  if(id %in% tko_pssea2$newID){
    # pheatmap function
    deseqs<- tko_deseq2[which(tko_deseq2$hgnc_symbol %in% c(gene1,gene2)),]
    deseq_tab<- deseqs[,c("sgrna","hgnc_symbol","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
    rawRPM<- tko_rawdata[deseq_tab$sgrna, which(grepl(cell, colnames(tko_rawdata)))]
    rawRPM2<- cbind(rawRPM[, which(grepl('Control', colnames(rawRPM)))], rawRPM[, which(grepl('Tcell', colnames(rawRPM)))])
    # for render a small table
    rownames(deseq_tab)<- NULL; colnames(deseq_tab)<- gsub('log2FoldChange','log2FC',colnames(deseq_tab));
    deseq_tab$Symbol<- data.frame(do.call('rbind', strsplit(deseq_tab$sgrna,'_',fixed=TRUE)))[,2]
    deseq_tab$sgrna<- data.frame(do.call('rbind', strsplit(deseq_tab$sgrna,'_',fixed=TRUE)))[,4]
    #deseq_tab<- deseq_tab[,c("Symbol","sgrna","baseMean","log2FC","lfcSE","stat","pvalue","padj")]
    deseq_tab<- deseq_tab[,c("Symbol","sgrna","log2FC","pvalue","padj")]
    
    anno_col<- data.frame(id=colnames(rawRPM2), grp=data.frame(do.call('rbind', strsplit(colnames(rawRPM2),'-',fixed=TRUE)))[,3])
    #anno_col$grp<- gsub('Control','Ctrl',anno_col$grp)
    anno_col<- data.frame(anno_col, row.names=1)
    anno_color = list(grp = c(Tcell="gold",Control="grey"))
    my_palette <- colorRampPalette(c("blue","white","red"))(n=29)
    pheat<- pheatmap(rawRPM2, col=my_palette, scale="row", cluster_row = F, cluster_col=F, 
                     treeheight_row = 0, #main="sgRNA profile log2(RPM + 1))", 
                     annotation_col=anno_col, annotation_colors = anno_color,
                     show_colnames = F, show_rownames = F, fontsize = 12, legend_labels = NULL,
                     gaps_row=c(as.numeric(table(deseqs$hgnc_symbol)[unique(deseqs$hgnc_symbol)[1]])))
    pheat<-  plot_grid(plot.new(), pheat[[4]],plot.new(), ncol=3, rel_widths=c(1,15,1))
    
    # barchat function charts
    pssea_out<- as.character(as.vector(tko_pssea2[which(tko_pssea2$newID == id),]))
    #pssea_out<- tko_pssea2[which(tko_pssea2$newID == id),]
    #print(pssea_out)
    pssea_tab<- data.frame(Genes=c(pssea_out[c(1:3)]),sgrna=c(pssea_out[c(6,4,5)]),
                           EnrichScore=c(pssea_out[c(7,9,11)]), pval=c(pssea_out[c(8,10,12)]))
    colnames(pssea_tab)<- c("Gene/s","# sgrna","Enrich Score","p-value")
    
    rankdata<- tko_deseq2[,c('sgrna','log2FoldChange','hgnc_symbol','NTC')]
    rankdata$Group<- ifelse(rankdata$hgnc_symbol %in% c(gene1, gene2), 'ParalogPair','REST')
    rankdata$labels<- ifelse(rankdata$NTC == 1, 'NTC',
                             ifelse(rankdata$hgnc_symbol == gene1, 'gene1',
                                    ifelse(rankdata$hgnc_symbol == gene2, 'gene2', 'Other')))
    rankdata2 = rankdata[which(rankdata$labels!='Other'),]
    rankdata2$size = ifelse(rankdata2$labels=='NTC','1','2')
    #table(rankdata$labels)
    ksplot<- ggplot(rankdata, aes(x=log2FoldChange, colour = Group)) + stat_ecdf() +
      scale_color_manual(values=c('#F8766D','lightgrey')) +
      theme_bw()
    
    cols <- c("gene1"="#2600D1","gene2"="#EE3F3F", 'NTC'='lightgrey')
    xmin = min(rankdata$log2FoldChange); xmax = max(rankdata$log2FoldChange);
    #rect= data.frame(x=c(xmin,xmin,xmax,xmax),y=c(0,1,0,1))
    rankplot = ggplot(data.frame(x=c(xmin,xmin,xmax,xmax),y=c(0,1,1,0))) + 
      geom_polygon(aes(x=x, y=y), fill="white", color="gray20") +
      geom_segment(aes(x=log2FoldChange,y=0, xend=log2FoldChange, yend=1, color = labels, size=size), data=rankdata2) +
      scale_color_manual(values = cols) +
      scale_size_manual(guide = "none", values = c(0.5, 2)) +
      scale_x_continuous(name="sgRNA ranked by log2FC", expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      ggtitle("Query paralog pair distributions") +
      theme(axis.line = element_blank(),
            plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y =element_blank(),
            axis.title.y = element_blank())
    plt<- plot_grid(rankplot, ksplot, ncol=1, rel_heights=c(1,3))
    plt2<- plot_grid(plot.new(), plt, plot.new(), ncol=3, rel_widths=c(1,15,1))
    
    return(list(g1_info, g2_info, deseq_tab, pheat, plt2, pssea_tab))
  }else{
    return(list(g1_info, g2_info)) # if genes not in the pSSEA results
  }
}


# Function for XGBoost result query
predict_visz_func<- function(gene1, gene2){
  #print(gene1)
  #print(gene2)
  id<- paste(min(gene1,gene2), max(gene1,gene2), sep='/')
  #print(id)
  #library(ggplot2)
  #library(ggrepel)
  #library(cowplot)
  
  if(id %in% pred$jointSymbol){
    rankdata<- pred[,c('jointSymbol','proba','ranking')]
    rankdata$cologrp<- factor(ifelse(rankdata$jointSymbol==id, 1, 0))
    
    # output 1
    options(ggrepel.max.overlaps = Inf)
    plt1<- ggplot(rankdata, aes(x = ranking, y = proba, color = cologrp)) +
      geom_point(aes(size=cologrp), alpha=0.6) +
      geom_text_repel(data=rankdata, aes(label=ifelse(cologrp==1, jointSymbol,'')),
                      box.padding = unit(2, "lines"), point.padding = unit(0.35, "lines"), 
                      segment.size  = 0.6, force_pull = 0,
                      segment.color = "grey50",
                      colour='orangered',size=4) + 
      scale_size_manual(values=c(0.5, 3)) +
      scale_color_manual(values=c("gray90",'#008bc6')) +
      labs(x = "Prediction rank", y = "Prediction score") +
      ggtitle('Ranking of XGBoost predicted probability') +
      theme_bw(base_size = 16) +
      theme(legend.position="none",
            plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"))+ 
      geom_vline(xintercept = round(dim(rankdata)[1] * 0.05), col = "gray50", linetype=3)
    
    # output 2
    featOut= data.frame(Feature=colnames(pred)[7:38], metric=as.numeric(pred[which(pred$jointSymbol==id), 7:38]))
    
    featDat= data.frame(feat=gsub('Zscore_','',colnames(pred)[41:72]), zscore=as.numeric(pred[which(pred$jointSymbol==id), 41:72]))
    shapDat= data.frame(feat=gsub('SHAP_','',colnames(pred)[73:104]), SHAP=as.numeric(pred[which(pred$jointSymbol==id), 73:104]))
    pltdata<- merge(featDat, shapDat, by='feat', all=F)
    pltdata$feat<- factor(pltdata$feat, levels=rev(pltdata$feat))
    pltdata$feat_cols<- ifelse(pltdata$zscore>1.2, "1",
                               ifelse(pltdata$zscore<=-1.2, "5",
                                      ifelse(pltdata$zscore<=1.2 & pltdata$zscore>0.4, "2",
                                             ifelse(pltdata$zscore<=0.4 & pltdata$zscore>-0.4, "3",
                                                    ifelse(pltdata$zscore<=-0.4 & pltdata$zscore>-1.2, "4", 0)))))
    
    p2<- ggplot(data=pltdata,  aes(x = SHAP, y = feat)) +
      geom_bar(stat = "identity", aes(fill = feat_cols)) +
      scale_fill_manual(breaks = c("1", "2", "3", "4", "5"),
                        values=c("#ff795e","#ffbfa4","#e4e3e3","#98d6ff","#3daaff"),
                        labels  = c("1", "2", "3", "4", "5"),
                        name = "Feature-Zscore") +
      xlab("SHAP value") + ggtitle("SHAP explaination for each feature") +
      theme_bw(base_size = 16) +
      theme(legend.position = 'none',
            plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y=element_blank(),
            axis.line = element_line(colour = "black"))
    
    p3<- ggplot(data.frame(y=1, x=c(1:5), grp=factor(c(1:5))), aes(x, y, fill = grp)) + 
      geom_bar(stat="identity",width=1) +
      scale_fill_manual(values=c("#3daaff","#98d6ff","#e4e3e3","#ffbfa4","#ff795e")) +
      scale_y_continuous(expand = c(0, 0)) +
      scale_x_continuous(expand = c(0, 0), breaks=c(),
                         sec.axis = sec_axis(trans = ~., breaks = c(1.5,2.5,3.5,4.5), labels = c('-1.2','-0.4','0.4',1.2))) +
      ggtitle("Feature−Zscore") +
      coord_flip() +
      theme_bw(base_size = 16) +
      theme(legend.position = 'none',
            plot.title = element_text(size=12),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x =element_blank(),
            axis.title =element_blank(),
            axis.line = element_blank())
    
    # output 3
    plt2<- plot_grid(p2, plot_grid(plot.new(), p3, plot.new(), ncol=1), plot.new(), ncol=3, rel_widths=c(8,1,0.5))
    plt3<- plot_grid(plt2, plot.new(), ncol=1, rel_heights=c(8,1))
    
    #print(featOut)
    return(list(plt1, plt3, featOut))
  }else{
    return(list(NA,NA,NA)) # if genes not in the results
  }
}

# Run the application
shinyApp(ui = ui, server = server)