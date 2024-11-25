library(MatchIt)
library(dplyr)
library(tableone)
library(reshape)
library(reshape2) 
library(ggpubr)

#import and arrange data --------
bray_dist <- read.table(paste(path.file, "bray_dist.txt",sep = ""),sep = "\t")
dist.m <- as.matrix(bray_dist)

for (i in dist_name) {
  print(i)
  dist.m = melt(get(i))
  print(class(dist.m))
  dist.m = dist.m %>%
    filter(as.character(Var1) != as.character(Var2)) %>%
    mutate_if(is.factor, as.character)

  index<- duplicated(paste0(pmax(dist.m$Var1, dist.m$Var2),pmin(dist.m$Var1, dist.m$Var2))) 
  dist.m <- dist.m[!index,]
  
  dist.m$pairs <- paste(dist.m$Var1,"_", dist.m$Var2, sep = "")
}

# 
data = dist.m

# 
randomID <- runif(nrow(data),0,1) 
r.data <- cbind(data, randomID) 
head(r.data)
r.data <- r.data[order(r.data$randomID),]
head(r.data)

# Propensity scores--------
# all related pairs
related <- data[which(r.data$over_related==0 | r.data$over_related==1),]
nrow(related)
# [1] 232903
attach(related)
m.related = matchit (over_related ~ avg_age + age_diff + gender_comb, method ="nearest", ratio =3)
m.related

f.related <-match.data(m.related, group = "all", distance = "distance", data = related, drop.unmatched = TRUE) 
write.table(f.related, paste(path.file, "matchPairs_related.txt", sep = ""),sep = "\t")

# first-relatives pairs
FirstRel <- r.data[-which(r.data$first_relative==2),]
FirstRel$first_relative <- as.factor(FirstRel$first_relative)
detach()
attach(FirstRel)
m.FirstRel = matchit (first_relative ~ avg_age + age_diff + gender_comb, method ="nearest", ratio =3)
f.FirstRel <-match.data(m.FirstRel, group = "all", distance = "distance", data = FirstRel, drop.unmatched = TRUE) #得到匹配后文件
head(f.FirstRel)
nrow(f.FirstRel)
write.table(f.FirstRel, paste(path.file, "matchPairs_FirstRel.txt",sep = "\t")

# parent-child pairs
Child <- r.data[-which(r.data$parent_offspring==2),]
nrow(Child)
detach()
attach(Child)
m.Child = matchit (parent_offspring ~ avg_age + age_diff + gender_comb, method ="nearest", ratio =3)
f.Child <-match.data(m.Child, group = "all", distance = "distance", data = Child, drop.unmatched = TRUE) #得到匹配后文件
head(f.Child)

write.table(f.Child, "matchPairs_Child.txt", sep = ""),sep = "\t")

# siblings
Siblings <- r.data[-which(r.data$siblings==2),]
Siblings$siblings <- as.factor(Siblings$siblings)
detach()
attach(Siblings)
m.Siblings = matchit (siblings ~ avg_age + age_diff + gender_comb, method ="nearest", ratio =3)
f.Siblings <-match.data(m.Siblings, group = "all", distance = "distance", data = Siblings, drop.unmatched = TRUE) 
write.table(f.Siblings, "ctrl+untreatedNPC_matchPairs_Siblings.txt",sep = "\t")

# spouses
Spouse <- r.data[-which(r.data$spouses==2),]
Spouse$spouse <- as.factor(Spouse$spouse)
detach()
attach(Spouse)
m.Spouse = matchit (spouse ~ avg_age + age_diff + gender_comb, method ="nearest", ratio =3)
f.Spouse <-match.data(m.Spouse, group = "all", distance = "distance", data = Spouse, drop.unmatched = TRUE) 
write.table(f.Spouse, "ctrl+untreatedNPC_matchPairs_Spouse.txt",sep = "\t")