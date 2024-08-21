options(java.parameters = "-Xmx25000m")
library(xlsx)
library(biomaRt)
library(progress)

ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
dataset_list<-datasets$dataset


############################################# USER INPUTS ##############################################################

filename<-"Zebrafish_SC_Cell type marker_database_052324_v7.xlsx"  #The filename of the xcel sheet where gene lists stored in separate columns
sheetname<-"Neurons_Anderson_mouse" #The sheet name inside the given excel file above
input_species<-"mmusculus" ##Type in the input species name, For e.g drerio for zebrafish, hsapiens for humans and mmusculus for mouse
species.to.convert<-"drerio"   ##Type in the orthologous species name, For e.g drerio for zebrafish, hsapiens for humans and mmusculus for mouse

##########################################################################################################################

data<- read.xlsx(file=filename, sheetName =sheetname )#Reading in the excel sheet with different list of genes in different columns
converted.data<-matrix(nrow = 2*(length(row.names(data))), ncol = length(colnames(data)), dimnames = list(NULL, colnames(data))) #creating matrix to store converted gene list

print(paste0("Potential Dataset that match the query for input species :", dataset_list[grep(input_species, dataset_list)]))
selected.input.dataset<-dataset_list[grep(input_species, dataset_list)][1]
print(paste("Selected input species dataset: ", selected.input.dataset))
print("Make sure this input species dataset is correct!!!")
input_dataset<- selected.input.dataset #If this is wrong, manually type in the exact dataset name that you found in potential Dataset that match the query for input species
ensembl.input <- useMart("ensembl", dataset = input_dataset)

# Creating a new progress bar to track time
pb <- progress_bar$new(
  format = " Finding orthologues genes [:bar] :percent eta: ",
  total = length(colnames(converted.data)), clear = FALSE, width= 100) 



  
#Defining gene attributes to be given inside the geBM function
gene_name<- paste0(species.to.convert,"_homolog_associated_gene_name")
gene_id<-  paste0(species.to.convert,"_homolog_ensembl_gene")

####checking if all the gene attributes are correct
attributes_list<-ensembl.input@attributes[,1]
sub_attributes<- attributes_list[grep(species.to.convert,attributes_list)]
is.element(gene_name,sub_attributes)
is.element(gene_id,sub_attributes)
  

##Starting the gene conversion using Biomart API
for (i in 1:length(colnames(converted.data))) {
    pb$tick()
    genes<- data[,i]
    genes<- unique(genes)
    genes<- genes[!is.na(genes)]
  
  
  ##searchFilters(mart = ensembl.input, pattern = "external_gene_name")
  
  ## Contacting biomart to identify orthologues
  ortho<- getBM(attributes = c("ensembl_gene_id","external_gene_name"  ,gene_name, gene_id), 
                filters ="external_gene_name",
                values = genes,
                mart = ensembl.input, uniqueRows = T)
  
  ortho_genes<-ortho[,3]
  ortho_genes<-unique(ortho_genes)
  ortho_genes<-ortho_genes[!is.na(ortho_genes)]
  ortho_genes<- ortho_genes[nchar(ortho_genes)>0]
  for (j in 1:length(ortho_genes)) {
    converted.data[j,i]<-ortho_genes[j]
  }
}


  
write.xlsx(converted.data, file =filename, append = T, sheetName = paste0(species.to.convert, " Converted"), showNA = F, row.names = F)



    

    
