# Trabalho de estatistica
# Abril de 2023
# Dados de https://www.kaggle.com/datasets/chaunguynnghunh/sepsis
# Julia Raspante Martins
# 
# ID: number to represent patient ID
# PRG: Plasma glucose
# PL: Blood Work Result-1 (mu U/ml)
# PR: Blood Pressure (mm Hg)
# SK: Serum potassium (mm)
# TS: Trasnferrin (mu U/ml)
# M11: Body mass index (weight in kg/(height in m)^2
# BD2: Beta-defensin 2 (mu U/ml)
# Age: patients age (years)
# Insurance: If a patient holds a valid insurance card
# Sepsis: Positive: if a patient in ICU will develop a sepsis , and Negative: otherwise

################################
library("dplyr")
library("ggplot2")
library("gridExtra")
library("grid")
library("PCAtools")
library("corrplot")
library("tydir")
library("vcd")

# importar arquivo de treinamento
# paitients_files_test <- read.csv("~/Paitients_Files_Test.csv")

patients <- read.csv("~/Paitients_Files_Train.csv", row.names=1)

# explorar o arquivo
names(patients) = c("PRG","PL","PR","SK","TS","M11","BD2","Age","Insurance","Sepsis")

# transformar Insurence e Sepsis em fator
patients$Insurance = as.factor(patients$Insurance)
patients$Insurance <- ifelse(patients$Insurance == 0, "no", "yes")

patients$Sepsis = as.factor(patients$Sepsis)

  
# removendo linhas com valores zero em qualquer coluna
patients <- subset(patients, !apply(patients == 0, 1, any))


patients_summary = as.data.frame(summary(patients))

# graficos descritivos

hist1 = function(dados, medida) {
  ggplot(dados, aes(x=dados[,medida])) +
    geom_histogram(bins = 15, fill="white", color="#919c9a") +
    geom_density(aes(y=..count..),fill="#779BE7", alpha=0.3, size=0.3, color="#779BE7") +
    labs(x = medida) +
    theme_minimal()
}


PRG = hist1(patients, "PRG")
PL = hist1(patients, "PL")
PR = hist1(patients, "PR")
SK = hist1(patients, "SK")
TS = hist1(patients, "TS")
M11 = hist1(patients, "M11")
BD2 = hist1(patients, "BD2")
Age = hist1(patients, "Age")

todos = grid.arrange(PRG,PL,PR,SK,TS,M11,BD2,Age, ncol=4)

barplot1 = ggplot(patients, aes(x=patients[,"Sepsis"],fill=Insurance)) +
  geom_bar(width = 0.5) +
  labs(x = "Sepsis") +
  theme_minimal() +
  scale_fill_manual(values = c("#779BE7","#F7717D"))

barplot1

boxplot2 = function(dados, medida) {
  ggplot(dados, color="#25282a") +
    geom_boxplot(aes(x=Sepsis, y=dados[,medida], fill=Sepsis),show.legend = FALSE) +
    labs(x = "", y=medida) +
    scale_fill_manual(values = c("#779BE7","#F7717D")) +
    theme_minimal()
}

PRG2 = boxplot2(patients, "PRG")
PL2 = boxplot2(patients, "PL")
PR2 = boxplot2(patients, "PR")
SK2 = boxplot2(patients, "SK")
TS2 = boxplot2(patients, "TS")
M112 = boxplot2(patients, "M11")
BD22 = boxplot2(patients, "BD2")
Age2 = boxplot2(patients, "Age")

todos2 = grid.arrange(PRG2,PL2,PR2,SK2,TS2,M112,BD22,Age2, ncol=4)


# teste de normalidade
shapiro.test(patients$PRG) # W = 0.8, p-value = 2e-15
shapiro.test(patients$PL) # W = 1, p-value = 5e-07
shapiro.test(patients$PR) # W = 1, p-value = 0.008
shapiro.test(patients$SK) # W = 1, p-value = 0.005
shapiro.test(patients$TS) # W = 0.8, p-value <2e-16  
shapiro.test(patients$M11) # W = 1, p-value = 0.04
shapiro.test(patients$BD2) # W = 0.9, p-value = 6e-14 
shapiro.test(patients$Age) # W = 0.9, p-value = 8e-15 

# nenhuma das as veriaveis testadas apresentam distribuicao normal

# separar grupo de sepsis positive e sepsis negative
patients_positive = patients |>
  filter(Sepsis == "Positive")

patients_negative = patients |>
  filter(Sepsis == "Negative")

# teste wilcoxon para amostras nao-pareadas
# H0: não ha diferenca na mediana entres grupos
# sepsis positive e sepsis negative
# H1: a mediana eh diferente para sepsis positive e sepsis negative
# alpha = 0.05
wilcox.test(patients_positive$PRG, patients_negative$PRG, paired = FALSE)
# W = 10130, p-value = 4e-06
wilcox.test(patients_positive$PL, patients_negative$PL, paired = FALSE)
# W = 12068, p-value = 2e-15
wilcox.test(patients_positive$PR, patients_negative$PR, paired = FALSE)
# W = 9851, p-value = 5e-05
wilcox.test(patients_positive$SK, patients_negative$SK, paired = FALSE)
# W = 10160, p-value = 4e-06
wilcox.test(patients_positive$TS, patients_negative$TS, paired = FALSE)
# W = 11812, p-value = 7e-14
wilcox.test(patients_positive$M11, patients_negative$M11, paired = FALSE)
# W = 10148, p-value = 5e-06
wilcox.test(patients_positive$BD2, patients_negative$BD2, paired = FALSE)
# W = 9129, p-value = 0.005
wilcox.test(patients_positive$Age, patients_negative$Age, paired = FALSE)
# W = 11612, p-value = 9e-13
# ao nivel de significancia de 5%, a mediana dos grupos sepsis positive
# e sepsis negative eh diferente em todas as variaveis testadas

# teste exato de Fisher para Insurance vs Sepsis
matriz_contingencia <- table(patients$Sepsis, patients$Insurance)
resultado_fisher = fisher.test(matriz_contingencia)
# p-value = 0.7 --> nao ha associacao entre sepsis e insurance
# 95 percent confidence interval:
#   0.63 2.09
# sample estimates:
#   odds ratio 
# 1.1

# excluir colunas de sepsis e insurance
patients_transpose = patients |>
  select(-c("Insurance", "Sepsis"))
# transpor dataframe para PCA
patients_transpose = as.data.frame(t(patients_transpose))

# criar dataframe de metadata para PCA
patients_metadata = data.frame(Sepsis = patients$Sepsis,
                      row.names = rownames(patients))

pca1 = pca(mat = patients_transpose,
           metadata = patients_metadata)

screeplot1 = screeplot(pca1, getComponents(pca1,components = 1:10),
                       title = "Scree plot")

biplot1 = biplot(pca1, x="PC1", y="PC2",
                 colby = "Sepsis",
                 colkey = c("#779BE7","#F7717D"),
                 legendPosition = "top")

# testar correlacao de Spearman
corr_spearman = cor(patients[,1:8], method = "spearman")

#corr_spearman_df <- reshape2::melt(corr_spearman)

spearman_corr <- read.delim("~/spearman_corr")

matriz_p <- cor.mtest(corr_spearman, conf.level = 0.95)$p

matriz_cor_p <- merge(spearman_corr, matriz_p, by = c("Var1", "Var2"))
names(matriz_cor_p) <- c("var1", "var2", "correlacao", "p_value")

spearman_corr_with_p <- read.delim("~/spearman_corr_with_p")

# Criando o heatmap com a função geom_tile()
corr_map = ggplot(spearman_corr_with_p, aes(var1, var2, fill = correlacao)) +
  geom_tile() +
  scale_fill_gradient2(low = "#779BE7",mid = "white", high = "#F7717D", midpoint = 0, limit = c(-1,1)) +
  coord_fixed() +
  theme_minimal() +
  geom_text(aes(label = round(correlacao, 2)), color = "black", size = 3) +
  labs(x=NULL, y=NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
corr_map

# graficos com as correlacoes encontradas
PRG_Age = ggplot(patients, aes(x=PRG, y=Age)) +
  geom_point(aes(color=Sepsis),show.legend = FALSE) +
  labs(x ="PRG", y="Age") +
  scale_color_manual(values = c("#779BE7","#F7717D")) +
  theme_minimal()

SK_M11 = ggplot(patients, aes(x=SK, y=M11)) +
  geom_point(aes(color=Sepsis),show.legend = FALSE) +
labs(x ="SK", y="M11") +
  scale_color_manual(values = c("#779BE7","#F7717D")) +
  theme_minimal()

TS_PL = ggplot(patients, aes(x=TS, y=PL)) +
  geom_point(aes(color=Sepsis),show.legend = FALSE) +
  labs(x ="TS", y="PL") +
  scale_color_manual(values = c("#779BE7","#F7717D")) +
  theme_minimal()

PR_BD2 = ggplot(patients, aes(x=PR, y=BD2)) +
  geom_point(aes(color=Sepsis),show.legend = FALSE) +
  labs(x ="PR", y="BD2") +
  scale_color_manual(values = c("#779BE7","#F7717D")) +
  theme_minimal()

PL_Age = ggplot(patients, aes(x=PL, y=Age)) +
  geom_point(aes(color=Sepsis),show.legend = FALSE) +
  labs(x ="PL", y="Age") +
  scale_color_manual(values = c("#779BE7","#F7717D")) +
  theme_minimal()

PL_SK = ggplot(patients, aes(x=PL, y=SK)) +
  geom_point(aes(color=Sepsis),show.legend = FALSE) +
  labs(x ="PL", y="SK") +
  scale_color_manual(values = c("#779BE7","#F7717D")) +
  theme_minimal()

PL_Age
PR_BD2
PL_SK

correlacoes = grid.arrange(PRG_Age, SK_M11, TS_PL,PL_Age,PL_SK,PR_BD2, ncol=3)



pdf("estatistica_figuras.pdf", width = 8, height = 4)
todos = grid.arrange(PRG,PL,PR,SK,TS,M11,BD2,Age, ncol=4)
barplot1
todos2 = grid.arrange(PRG2,PL2,PR2,SK2,TS2,M112,BD22,Age2, ncol=4)
corr_map
correlacoes = grid.arrange(PRG_Age, SK_M11, TS_PL,PL_Age,PL_SK,PR_BD2, ncol=3)
dev.off()
