


#_______________________________________________________________________________________________________
#Her er primerne slik de står i labheftet vårt nanopore sekvensering.  



#[1] Modified Mangala-F1 primer: 5’ TTTCTGTTGGTGCTGATATTGC TCCTACGGGAGGCAGCAG 3’



#[1] Modified 16S UR primer: 5’ ACTTGCCTGTCGCTCTATCTTC- CGGTTACCTTGTTACGACTT 3’

--------------------------------------------------------------------------------------------------------
  
#Vedlagt er et utvalg av 16S sekvenser fra SILVA-databasen, i form av en fasta-fil. 
#Jeg har en match mot en sekvens dersom den ene primeren matcher ett sted, og den andre primeren matcher et annet sted. Merk at begge primere er angitt i 5’-3’ retning, det betyr vel at den ene må revers-komplementeres før man matcher. 
#Jeg er interessert i antall sekvenser som matcher, samt deres taksonomi. Det er angitt i Header-linjen til hver sekvens, adskilt med semikolon for hvert nivå. Dette følger stort sett mønsteret:
#Superkingdom;phylum;class;order;family;genus;species
#Deretter jeg ønsker å ha listet genus hos alle sekvenser med primer-match.


#Den vedlagte fasta-fila er bare et utplukk på 1000 sekvenser fra den store SILVA samlingen. Når koden foreligger, last ned fila SILVA_138_SSURef_tax_silva.fasta.gz fra https://www.arb-silva.de/no_cache/download/archive/release_138/Exports/
  #Archive - ARB & SILVA




library(tidyverse)
library(dplyr)


lines <- readLines("SILVA_138_SSURef_tax_silva_sub1000.fasta")

header.idx <- which(str_detect(lines, ">")) #detektere linjer som inne holder ">"



fasta.table <- tibble(Header = str_remove(lines[header.idx], "^>"),
                      Sequence = rep("", length(header.idx)))

print(fasta.table)

header.idx <- c(header.idx, length(lines)+1)
for(i in 1:(length(header.idx)-1)){
  chunk <- lines[(header.idx[i]+1):(header.idx[i+1]-1)]
  fasta.table$Sequence[i] <- str_c(chunk, collapse = "")  
}





#________________________________________________________________________________________________________

library(stringr)
library(tidyverse)
library(tidyr)


genus <- fasta.table$Header

genus.table <- c(genus)
b <- list(strsplit(genus.table, ";"))

df <- data.frame(genus = (fasta.table$Header)) 



library(splitstackshape)

genusfinal.table <- cSplit(data.frame(genus), 'genus', ";")

genusfinal.table[,`:=`(genus_1 = NULL, genus_2= NULL, genus_3= NULL, 
                       genus_4= NULL, genus_5= NULL, genus_7= NULL)] #genusfinal.table Printer oppdaterte genusfinal.table



RNA <- c(fasta.table$Sequence)
riktig <- gsub("U", "T", RNA)


#Konverterer til DNA.
sekvenser.table <- riktig




#________________________________________________________________________________________________________

#Oppdaterte tabellen


final.table.ext <- cbind(genusfinal.table, sekvenser.table)   #legger sammen sekvenser og genus tabellen sammen       





final.table.ext%>% slice(400:1000) %>% 
  separate(sekvenser.table, into = c("Sekvens", NA), sep = 10)






#________________________________________________________________________________________________________


#Legger til primeren slik som det står i labheftet vårt nanopore sekvensering.Deretter setter final.table.ext og primer sammen.


# Prøver å mathche primer sekvensen med den oppdaterte tabellen.



library(stringr)



#primer <- data.frame(Modified_Mangala_F1_primer = c("TTTCTGTTGGTGCTGATATTGC TCCTACGGGAGGCAGCAG "), 
#Modified_16S_UR_primer  = c("ACTTGCCTGTCGCTCTATCTTC- CGGTTACCTTGTTACGACTT "), 
#stringsAsFactors = FALSE)

total.table.ext <- cbind(final.table.ext, primer)

#5´-3
res_1 = final.table.ext %>%
  filter(str_detect(final.table.ext$sekvenser.table, "TCCTACGGGAGGCAGCAG"))






res_1 %>% slice(800:901) %>% 
  separate(sekvenser.table, into = c("Sekvens", NA), sep = 10) 





#primer_2 3´-5´ og KOMPLEMENTÆR!!!!
res_2 = final.table.ext %>%
  filter(str_detect(final.table.ext$sekvenser.table, "AAGTCGTAACAAGGTAACCG"))




res_2 %>% slice(20:40) %>% 
  separate(sekvenser.table, into = c("Sekvens", NA), sep = 10) 

#begge to
res_3 = final.table.ext %>%
  filter(str_detect(final.table.ext$sekvenser.table, "TCCTACGGGAGGCAGCAG")&(str_detect(final.table.ext$sekvenser.table, "AAGTCGTAACAAGGTAACCG")))



res_3 %>% slice(1:10) %>% 
  separate(sekvenser.table, into = c("Sekvens", NA), sep = 10) 

#husk å kommentar at du har 39 matcher
#

#Du kan teste om denne koden fungerer. Da vil du først velge ut kun de fem første radene med Slice(1:5). #Deretter vil du dele kolonnen KOLONNENAVN inn i to deler. Du deler etter 20 tegn. Den første delen #havner i en ny kolonne kalt Sekvens, mens resten "forkastes". Det er da viktig å fortelle leseren at du #slik bare viser en liten del av res_3, og at den egentlig inneholder flere rader, og at DNA-sekvensen #egentlig er mye lengre.

#slice er fra dplyr-pakka
#seperate er fra tidyr-pakka


#_______________________________________________________________________________________________________

#Den første primer må være i den retningen oppgitt i protokollen, men den andre må være i 3´- 5´ og komplementær.

#Når jeg bruker bare første får jeg 901 match, deretter bruker jeg primer nr2 og da får jeg 40. 
#Disse to primere skal festes til og amplifisere et området da trenger du å bruke filter som har begge sekvenser.
#Dette fører til at jeg får 39 hits. Dvs at disse primere kan brukes på 39 av 1000 genomer.

#______________________________________________________________________________________________________
# Variasjon i GC innhold.




library(seqinr)
DNA <- read.fasta(file = "SILVA_138_SSURef_tax_silva_sub1000.fasta")
DNAsek <- DNA[[1]]

GC(DNAsek)






#I molekylærbiologi  er GC innhold prosentandelen av nitrogenholdige baser i et DNA- eller RNA molekyl som er enten guanin eller cytosin. GC innheoldet i 1000 første i Silva dna sekvensen er omtrent 65,5 %, er det sannsynligvis lokal variasjon i GC innhold i genomet. Det er noen regioner i genomsekvensen kan ha GC innhold ganske mye høyere enn 65,5%, mens noen regioner i genomssekvensen kan ha GC innhold som er ganske mye lavere enn 
#65,5% .



starts <- seq(1, length(DNAsek)-1000, by = 1000)
starts

n <- length(starts)    # Find the length of the vector "starts"
for (i in 1:n) {
  chunk <- DNAsek[starts[i]:(starts[i]+999)]
  chunkGC <- GC(chunk)
  print (chunkGC)
}




## Graf av GC innhold


starts <- seq(1, length(DNAsek)-10, by = 1000)
n <- length(starts)    #  Skal finne lengde av en vektor "start"
chunkGCs <- numeric(n) # Make a vector of the same length as vector "starts", but just containing zeroes
for (i in 1:n) {
  chunk <- DNAsek[starts[i]:(starts[i]+999)]
  chunkGC <- GC(chunk)
  print(chunkGC)
  chunkGCs[i] <- chunkGC
  
}

plot(starts,chunkGCs,type="b",xlab="Nucleotide start position",ylab="GC content")


#Lokale svigninger i GC innhold i genomssekvensen kan gi forskjelllige interessant informasjon.
#Man kan for eksempel avsløre tilfeller av horisontal overføring eller forstyrrelser i mutasjoner.




GCpros <- double()
n <- length(DNA[[1]])
for (i in 1:(n - 999)) GCpros[i] <- GC(DNA[[1]][i:(i+999)])
plot(GCpros,type="l")



#Hvis en del av DNA har beveget seg ved horisontal overføring fra genomet til lav GC innhold til med høyt GC innhold, kan chunk  av horisontalt overført oppdages som en region med uvanlig lavt GC innhold i høy GC mottakerens genom.

#Generelt er det slik at lav GC innhold i et ellers høy GC innholdsgenom  kan også oppstå av forstyrrelser i mutasjon i den regionen av genomet, for eksempel hvis mutasjoner fra G-C til T-A er vanligere av en eller annen grunn i den regionen av genomet enn i resten av genomet.
#________________________________________________________________________________________________________

#Nå ønsker jeg å finne  hvilke phyla andelen matcher et større/mindre forventet.


library(stringr)
library(tidyverse)
library(tidyr)


genus <- fasta.table$Header

genus.table <- c(genus)
b <- list(strsplit(genus.table, ";"))

df <- data.frame(genus = (fasta.table$Header)) 



library(splitstackshape)

genus.table <- cSplit(data.frame(genus), 'genus', ";")

genus.table[,`:=`(genus_1 = NULL, genus_3= NULL, 
                  genus_4= NULL, genus_5= NULL,genus_6= NULL, genus_7= NULL)] #genusfinal.table Printer oppdaterte genusfinal.table



RNA <- c(fasta.table$Sequence)
riktig <- gsub("U", "T", RNA)

# Loop through each row of 'dt' to replace Synonyms with word using sapply




sekvenser.table <- riktig

phylum.table.ext <- cbind(genus.table, sekvenser.table)   #legger sammen sekvenser og genus tabellen sammen       





#Her har jeg laget en tabelloversikt over alle phylum og sekvensene.


phylum.table.ext %>% slice(1:10) %>% 
  separate(sekvenser.table, into = c("Sekvens", NA), sep = 10) 
#________________________________________________________________________________________________________



phyl_res_1 = phylum.table.ext %>% 
  filter(str_detect(phylum.table.ext$sekvenser.table, "TCCTACGGGAGGCAGCAG"))


table(phyl_res_1)


#________________________________________________________________________________________________________




library(ggplot2)

ggplot(data = phyl_res_1,
       mapping = aes(x = genus_2)) + 
  geom_bar() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(y= "Antall sekvenser", x = "Phylum") 




#________________________________________________________________________________________________________
#Når jeg prøver å matche med primer Modified_Mangala_F1_primer, får 900 matcher ser vi at det er Phylum Firmicutes og Proteobacteria som dominerer.






#________________________________________________________________________________________________________


#primer_2 3´-5´ og KOMPLEMENTÆR!!!!
phyl_res_2 = phylum.table.ext %>%
  filter(str_detect(phylum.table.ext$sekvenser.table, "AAGTCGTAACAAGGTAACCG"))


#________________________________________________________________________________________________________


#finner oversikt over hvor hvilke phylum som ble matchet, og det var totalt  40  som var matchet av 1000 sekvenser. Utifra grafen nede ser vi at Proteobakterier dominerer.


ggplot(data = phyl_res_2,
       mapping = aes(x = genus_2)) + 
  geom_bar() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(y= "Antall sekvenser", x = "Phylum")  




#_____________________________________________________________________________________________________




#Når jeg prøver å matche med primer Modified 16S UR primer , får 40 matcher ser vi at det er Phylum Proteobacteria som dominerer.
#________________________________________________________________________________________________________

#begge to
phyl_res_3 = phylum.table.ext %>%
  filter(str_detect(phylum.table.ext$sekvenser.table, "TCCTACGGGAGGCAGCAG")&(str_detect(final.table.ext$sekvenser.table, "AAGTCGTAACAAGGTAACCG")))

#________________________________________________________________________________________________________




#Finner oversikt over hvor hvilke phylum som ble matchet, og det var totalt  39 som var matchet.


table(phyl_res_3)
#________________________________________________________________________________________________________

ggplot(data = phyl_res_3,
       mapping = aes(x = genus_2)) + 
  geom_bar() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(y= "Antall sekvenser", x = "Phylum")  





#_____________________________________________________________________________________________________




#Fra table(phyl_res_1) får jeg oversikt over hvor mange forskjellige phylum ble matchet.

df <- data.frame(
  Acidobacteriota = c(12,12),
  Actinobacteriota = c(135,138),
  Bacteroidota = c(90,102),
  Bdellovibrionota = c(1,1),
  Campilobacterota = c(8,8),
  Chloroflexi = c(5,10),
  Cloacimonadota = c(1,1),
  Cyanobacteria  = c(26,28),
  Dadabacteria = c(1,1),
  Deinococcota = c(1,1),
  Dependentiae = c(1,1),
  Desulfobacterota = c(13,13),
  Firmicutes = c(302,316),
  Fusobacteriota = c(1,1),
  Gemmatimonadota  = c(1,1),
  SAR406_clade = c(1,1),
  Methylomirabilota = c(1,1),
  Myxococcota  = c(3,3),
  NB1_j = c(1,1),
  Nitrospinota = c(1,1),
  Nitrospirota = c(3,3),
  Patescibacteria = c(2,4),
  Proteobacteria = c(282,288),
  SAR324_clade_Marine_group_B = c(1,1),
  Spirochaetota = c(5,5),
  Sva0485 = c(1,1),
  Synergistota = c(1,2),
  row.names = c('match',' ikke Match '))

df

fisher.test(df, simulate.p.value = TRUE)




#Nullhypotesen er  sannsynligheten for suksess i denne eksempelen så handler det om  Match og ikke match. Men problemet er at denne felles sannsynligheten er ukjent selv om nullhypotesen er sann og vi trenger for å beregne p-verdien. I oppgaven har jeg ikke så mye tall da ønsker jeg å  forholde meg til Fisher eksakt test,  ved å late som om det totale match er bestem før.

#Jeg får en p-verdi = 1 som betyr at jeg ikke klarte å avvise H0, dvs , det er ingen postivt sammenheng mellom variablene dine. Fisher eksakt test anbefales for små utvalgstørrelser og prøvestørrelsern som er relativt stor noen områder.
______________________________________________________________________________________________________

#table(phyl_res_2$genus_2)

ef <- data.frame(
  Actinobacteriota = c(1,138),
  Bacteroidota = c(2,102),
  Campilobacterota = c(4,8),
  Chloroflexi = c(5,10),
  Crenarchaeota = c(1,11),
  Cyanobacteria  = c(1,28),
  Firmicutes = c(1,316),
  Patescibacteria = c(2,4),
  Proteobacteria = c(30,288),
  row.names = c('match',' ikke Match '))

ef

fisher.test(ef, simulate.p.value = TRUE)
```



#table(phyl_res_3$genus_2)



ff <- data.frame(
  Actinobacteriota = c(1,138),
  Bacteroidota = c(2,102),
  Campilobacterota = c(4,8),
  Cyanobacteria  = c(1,28),
  Firmicutes = c(1,316),
  Proteobacteria = c(30,288),
  row.names = c('match',' ikke Match '))

ff

fisher.test(df, simulate.p.value = TRUE)
```


#________________________________________________________________________________________________________





## Konklusjon

#Jeg har funnet ut at primerne "Modified Mangala- F1 primer" og "Modified 16s UR primer" matcher med de sekvensene. Begge primere ble angitt i 5'-3' retning, "Modified 16s UR primer" ble revert komplementert.I denne oppgaven fikk jeg bare brukt et utplukk på 10000 sekvenser .fasta fil  fra den store SILVA samlingen. I følge mine analyser ble flertallet ble annotert til phylum Proteobacteria og Firmicutes. Proteobakterier er en gram negativ bakterie og Firmicutes er en undergruppe av gram positive bakterier.
#Fra resultater phyl_res_1 , phyl_res 2 og  Phyl_res_3 får jeg en p-verdi = 1 og klarte ikke å å avvise H0, dvs , det er ingen postivt sammenheng mellom variablene dine. 
#Den store SILVA samlingen består av drøyt 2 millioner sekvenser, som dere ser under har jeg prøvd å kjøre dette og har ventet i over 2 dager og den ble aldri ferdig kjørt. Jeg vet ikke helt om det er noe feil med maskinen eller om at jeg har en treig maskin. Har prøvd å bruke readFasta og readseqir, men det har ikke funket.
#Dette ønsker jeg å fortsette på min masteroppgaven.


#library(microseq)

#Silva_2mill <- readFasta("SILVA_138_SSURef_tax_silva.fasta")

#Silva.idx <- which(str_detect(Silva_2mill, ">")) #detektere linjer som inne holder ">"


#Silva_2mill.table <- tibble(Header = str_remove(Silva_2mill[Silva.idx], "^>"),
#Sequence = rep("", length(Silva.idx)))

#header.idx <- c(Silva.idx, length(Silva_2mill)+1)
#for(i in 1:(length(Silva.idx)-1))
#chunk <- Silva_2mill[(Silva.idx[i]+1):(Silva.idx[i+1]-1)]
#Silva_2mill.table$Sequence[i] <- str_c(chunk, collapse = "")  











