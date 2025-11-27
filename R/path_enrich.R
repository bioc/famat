#' @export path_enrich
#############################FUNCTIONS###########################
##perform pathway enrichment analysis on user's gene list
gene_path<-function(gene, db){

  entrez_ids <- clusterProfiler::bitr(gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID

  if(db == "KEGG"){
    res <- clusterProfiler::enrichKEGG(
      gene         = entrez_ids,
      organism     = 'hsa',     # Homo sapiens
      pvalueCutoff = 1,         # correspond à significant = FALSE
      qvalueCutoff = 1
    )
  }
  else if(db == "REAC"){
    res <- ReactomePA::enrichPathway(
      gene         = entrez_ids,
      organism     = "human",
      pvalueCutoff = 1,  # pour récupérer toutes les voies
      qvalueCutoff = 1,
      readable     = TRUE # convertit les Entrez IDs en symboles
    )
  }
  else if(db == "WP"){
    res <- clusterProfiler::enrichWP(
      gene         = entrez_ids,     # Symboles
      organism     = "Homo sapiens",
      pvalueCutoff = 1,
      qvalueCutoff = 1
    )
  }
  # Extraire id et description
  df_res <- as.data.frame(res)[, c("ID","Description")]
  colnames(df_res) <- c("id","name")
  df_res_filter<- subset(df_res, !is.na(df_res$id))

  return(df_res_filter)
}

##perform pathway enrichment analysis on user's metabolites list
meta_path<-function(meta, db){
    ##convert kegg ids in pubchem ids
    meta_pubchem<-mapper[which(mapper[, 1] %in% meta), 2]

    ##perform pathway enrichment analysis on user's metabolites list
    ##extract informations about pathways from MPINet results
    #pss<-getPSS(meta_pubchem, plot=FALSE)#removed MPINet::
    if(db == "KEGG"){
        medium<-identifypathway(meta_pubchem, "KEGG")#removed MPINet::
        resm<-vapply(medium, function(x){x[[2]][[1]]}, character(1))
    }else if(db == "REAC"){
        medium<-identifypathway(meta_pubchem, "Reactome")#removed MPINet::
        resm<-vapply(medium, function(x){x[[1]][[1]]}, character(1))
    }else if (db == "WP"){
        medium<-identifypathway(meta_pubchem, "Wikipathways")#removed MPINet::
        resm<-vapply(medium, function(x){x[[1]][[1]]}, character(1))
    }

    resm=as.data.frame(resm);names(resm)="name"
    return(resm)
}

##put pathway enrichment analysis results together and add pathways ids
path_enrich<-function(source, metabo, genes){
    ##perform pathway enrichment analysis on user's gene list
    resgene<-gene_path(genes, source)
    ##perform pathway enrichment analysis on user's metabolites list
    resmeta<-meta_path(metabo, source)
    ##add pathways ids
    resmeta<-cbind(resmeta, id=rep(NA, nrow(resmeta)))
    if(source == "KEGG"){
        prev_db<-KEGGREST::keggList("pathway", "hsa")
        idskegg=sub('^hsa:?(.*)$','\\1',names(prev_db))##Modif MC 11/04/23 => change in KEGGrest-Kegg pathways names?
        #idskegg=stringr::str_sub(names(prev_db), 9, nchar(names(prev_db)))
        idskegg=paste("hsa:",idskegg,sep="")
        keggdb<-as.list(idskegg)
        nameskegg<-unlist(stringr::str_split(unname(prev_db),
                                            " - Homo sapiens (.)human(.)"))
        names(keggdb)<-nameskegg[!(nameskegg %in% "")]

        resgene$id <- paste0("hsa:", sub("^hsa:?", "", resgene$id))
        resmeta[,2]=apply(resmeta, 1, function(x){
            id=keggdb[[x[1]]]
            if(!is.null(id)){id}
            else{NA}
        })
    }
    else if(source == "REAC"){
        readb<-unlist(as.list(reactome.db::reactomePATHID2NAME))
        readb<-readb[which(stringr::str_sub(readb, 1, 14)%in%"Homo sapiens: ")]
        pathids=names(readb);readb<-stringr::str_sub(readb, 15, nchar(readb))
        names(readb)=pathids
        resmeta[,2]=apply(resmeta, 1, function(x){
            id=names(readb[which(readb %in% x[1])])
            if(length(id)>0){id}
            else{NA}
        })
    }
    else if (source == "WP"){
        wpdb<-rWikiPathways::listPathways("Homo sapiens")
        wpdb<-wpdb[, c(3, 1)]
        ids=apply(resmeta, 1, function(x){
            pathname=stringr::str_to_upper(x[1])
            id<-wpdb[which(pathname == stringr::str_to_upper(wpdb$name)), 2]
            if(length(id)>0){id[1]}
            else{NA}
        })
        resmeta[,2]=unname(unlist(ids))
        resgene$id<-paste("WP", resgene$id, sep="")
    }
    resgene_filter <- resgene[!is.na(resgene$id), ]
    resmeta_filter <- resmeta[!is.na(resmeta$id), ]

    return(list(resmeta_filter, resgene_filter, genes, metabo, resmeta_filter))
}
