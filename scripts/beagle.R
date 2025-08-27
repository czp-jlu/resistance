#######################################
# HaploSelekt: Cross-validation study #
# Carola Zenke-Philippi               #
#######################################
# Imputing marker data with BEAGLE    #
#######################################

library(SelectionTools)
st.input.dir <- "input"
st.output.dir <- "output"

marker <- read.table("input/haploselekt-marker_preselected.mpo", header=T)
# preselected markers according to thresholds and place in input folder

# Recode marker data
# 1: A
# 2: C
# 3: G
# 4: T

marker[which(marker=="1/1", arr.ind=TRUE)] <- "A/A"
marker[which(marker=="2/2", arr.ind=TRUE)] <- "C/C"
marker[which(marker=="3/3", arr.ind=TRUE)] <- "G/G"
marker[which(marker=="4/4", arr.ind=TRUE)] <- "T/T"
marker[which(marker=="-1/-1", arr.ind=TRUE)] <- "0/0"

# order map and marker data by map position
map <- read.table("input/haploselekt-map.txt", header=T) # map must be in input folder
map <- map[order(map$pos),]
map <- map[order(map$chrom),]

# only markers in the map which are present in the data set
map <- map[map$name %in% rownames(marker),]

# are there duplicate map positions?
map[duplicated(map$pos),]
# yes
map$pos[duplicated(map$pos)] <- map$pos[duplicated(map$pos)] + 0.000001
# are there still duplicate map positions?
map[duplicated(map$pos),]
# yes
map$pos[duplicated(map$pos)] <- map$pos[duplicated(map$pos)] + 0.000001
# still?
map[duplicated(map$pos),]
# no

write.table(map, "input/haploselekt-map_sorted.txt",
            quote=F, row.names=F, col.names=F)

marker[1:5,1:5]
m <- marker[map$name,]
gen <- colnames(m)
m <- t(m)

# Read data into SelectionTools

st.read.marker.data ( filename="haploselekt-marker_preselected.mpo",
                     format="m")
st.read.map("haploselekt-map_sorted.txt", # Read in linkage map, here in Mbp
            m.stretch=1, # Multiplier for map positions
            format="mcp", # Marker Chromosome Position / "cpms"
            skip=1) # Skip the first line

###############################################################################
# Generate ped file for plink
###############################################################################
m <- cbind(gen,m)
colnames(m) <- NULL
cat("",file="haploselekt-plink.ped")
for (i in (1:nrow(m)))
write( gsub('/',' ', paste(m[i,],collapse=" " )),
      "haploselekt-plink.ped", append=TRUE)

ped <- read.table("haploselekt-plink.ped")

###############################################################################
# Generate map file for plink
###############################################################################
mp <- map[,c(2,1,3,3)]
mp[,3] <- mp[,3]/100 # Genetic map in M
mp[,4] <- mp[,4]*1000000 # Base pair position 1 Mbp - 1 cM
unlink("haploselekt-plink.map")
options(scipen=100) # Avoid scientific notation in the map (1e8)
write.table(mp,"haploselekt-plink.map",quote=F,row.names=F,col.names=F)
cmd <- paste ("tools/plink --file haploselekt-plink",
              "--no-fid --no-parents --no-sex --no-pheno --alleleACGT",
              "--recode vcf --out haploselekt-beagle")

# For Windows
# cmd <- paste ("tools//plink.exe --file v13-u03-3-plink",
# "--no-fid --no-parents --no-sex --no-pheno --alleleACGT",
# "--recode vcf --out v13-u03-3-beagle")
system(cmd)

###############################################################################
# Phasing with beagle
###############################################################################
cmd <- paste ("java -Xmx8g -jar tools/beagle.05May22.33a.jar",
              "nthreads=4 ap=true gp=true",
              "gt=haploselekt-beagle.vcf",
              "out=haploselekt-beagle-out")
# For Windows
# cmd <- paste ("java -Xmx8g -jar tools//beagle.08Nov19.3ec.jar",
# "nthreads=4 ap=true gp=true",
# "gt=v13-u03-3-beagle.vcf",
# "out=v13-u03-3-beagle-out")
b.captured <- system(cmd,intern=TRUE)
n <- length(b.captured)
b.captured[c(1:15,(n-5):n)]

system("gunzip -f haploselekt-beagle-out.vcf.gz")
# For Windows
# system("tools//gzip -d -f v13-u03-3-beagle-out.vcf.gz")
bd <- read.table("haploselekt-beagle-out.vcf",stringsAsFactors=F)
bd[1:17,1:18]

###############################################################################
# Format the Beagle output for SelectionTools
###############################################################################
m.matrix <- scan("haploselekt-beagle-out.vcf",character(),skip=11,nlines=1 )
m.matrix <- as.character(m.matrix[10:length(m.matrix)])
m.matrix <- matrix(unlist(strsplit(m.matrix,"_")),ncol=2,byrow=TRUE)[,1]
m.matrix <- gsub("\\.1","",m.matrix)
m.matrix <-paste(m.matrix,collapse=" ")
for (i in 1:nrow(bd))
{
m <- paste(as.character(bd[i,10:ncol(bd)]),collapse=" ")
m <- gsub ( "0", bd[i,4], m)
m <- gsub ( "1", bd[i,5], m)
m <- gsub ( "\\|", "/", m)
m <- gsub ( "\\.", "-", m)
m <- paste (as.character(bd[i,3]),m)
m.matrix <- rbind(m.matrix,m)
}
write.table (m.matrix,"haploselekt-beagle-out.mpo", quote=FALSE,
             row.names=FALSE,col.names=FALSE)

dim(m.matrix)
