# library(dplyr)
# library(ggplot2)
# 
# 
# # Définir le chemin du dossier
# folder_path <- "/Users/wenjiehuang/Desktop/LyonLumiere2/M2/StageM2/Data/DataComp"
# # Lister tous les fichiers CSV
# csv_files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)
# 
# # Initialiser un dataframe vide
# comp_data <- data.frame()
# 
# # Définir le taux d'échantillonnage d'origine et le taux d'échantillonnage cible
# original_sampling_rate <- 0.002
# target_sampling_rate <- 0.05
# 
# # Calculer le nombre de lignes à conserver
# rows_to_keep <- round(target_sampling_rate / original_sampling_rate)
# 
# # Parcourir chaque fichier CSV
# for (file in csv_files) {
#   # Extraire le numéro de sujet et le scénario expérimental du nom de fichier
#   file_name <- basename(file)
#   subject_number <- strsplit(file_name, "_")[[1]][1]
#   scenario <- strsplit(file_name, "_")[[1]][2]
#   
#   # Lire le fichier CSV
#   csv_data <- read.csv(file, header = TRUE, sep = ";", dec = ",")
#   
#   # Ajouter des colonnes subject_number et session
#   csv_data$Subject <- subject_number
#   csv_data$Session <- scenario
#   
#   # Échantillonner les données
#   sampled_data <- csv_data[seq(1, nrow(csv_data), by = rows_to_keep), ]
#   
#   # S'assurer que TypeCond et NbTrial sont des types appropriés
#   sampled_data$TypeCond <- as.character(sampled_data$TypeCond)
#   sampled_data$NbTrial <- as.numeric(sampled_data$NbTrial)
#   
#   # Ajouter TrialOrder en fonction du scénario expérimental et du numéro de sujet
#   sampled_data <- sampled_data %>%
#     group_by(Subject, Session) %>%
#     arrange(Time) %>%
#     mutate(
#       TrialOrder = 1 + cumsum(if_else(TypeCond != lag(TypeCond, default = first(TypeCond)) |
#                                         NbTrial != lag(NbTrial, default = first(NbTrial)), 1, 0))
#     )%>% ungroup()
#   # mutate(
#   #   TrialOrder = cumsum(c(TRUE, lag(TypeCond, default = first(TypeCond)) != TypeCond |
#   #                           lag(NbTrial, default = first(NbTrial)) != NbTrial))
#   # )
#   
#   # Ajouter les données au dataframe total
#   comp_data <- bind_rows(comp_data, sampled_data)
# }
# 
# # Remplacer toutes les valeurs de 2 par 0,5 dans BinaryChoice
# comp_data <- comp_data %>%
#   mutate(BinaryChoice = ifelse(BinaryChoice == 2, 0.5, BinaryChoice))
# 
# 
# # Effectuer d'abord la conversion de mappage
# comp_data_modified <- comp_data %>%
#   mutate(
#     BinaryChoice = 2 * as.numeric(BinaryChoice) - 1,
#     JoystickYAxis = 2 * JoystickYAxis - 1
#   )
# 
# View(comp_data_modified)
# 
# 
# comp_data_modified <- comp_data_modified %>%
#   group_by(Subject, Session, TrialOrder) %>%
#   mutate(
#     # Obtenir la dernière valeur de BinaryChoice pour chaque bloc de données
#     last_binary_choice = last(BinaryChoice),
#     last_RT = last(RT)
#   ) %>%
#   filter(last_binary_choice != 0)
# 
# 
# ####################### confirmer la corrélation entre joystick et choix binaire ####################
# ########################### comparer la corrélation entre les sessions ###############################
# 
# # Tracer la courbe de régression linéaire
# ggplot(comp_data_modified, aes(x = JoystickYAxis, y = as.numeric(last_binary_choice), color = Session)) +
#   geom_smooth(method = "lm", se = TRUE) + # Tracer la courbe de régression linéaire et l'intervalle de confiance
#   labs(title = "Relation entre JoystickYAxis et BinaryChoice par session",
#        x = "Axe Y du Joystick",
#        y = "Choix Binaire") +
#   theme_minimal() +
#   theme(legend.position = "bottom")
# 
# 
# # Ajuster le modèle linéaire combiné avec le terme d'interaction
# combined_lm <- lm(as.numeric(last_binary_choice) ~ JoystickYAxis * Session, data = comp_data_modified)
# # Résumé du modèle combiné
# summary(combined_lm)
# # ANOVA pour tester la signification de l'interaction
# anova(combined_lm)
# 
# ######################## coefficient de corrélation de Kendall ######################
# 
# # Définir la fonction pour calculer tau de Kendall et la p-value
# calculate_kendall_tau <- function(data) {
#   result <- cor.test(data$JoystickYAxis, as.numeric(data$last_binary_choice), method = "kendall")
#   return(data.frame(kendall_tau = result$estimate, p_value = result$p.value))
# }
# 
# # Calculer tau de Kendall et la p-value pour chaque session
# session_stats <- comp_data_modified %>%
#   group_by(Session) %>%
#   do(calculate_kendall_tau(.))
# 
# print(session_stats)
# 
# ######################### Modèle à effets mixtes sur RT~BinaryChoice ############
# 
# library(lme4)
# 
# # Modèle à effets mixtes : relation entre RT et BinaryChoice, considérant le sujet comme un effet aléatoire
# rt_mixed_model <- lmer(last_RT ~ last_binary_choice * Session + (1|Subject), data = comp_data_modified)
# summary(rt_mixed_model)
# 
# 
# ### RT par BinaryChoice à travers les sessions
# 
# summary_data <- comp_data_modified %>%
#   group_by(Session, last_binary_choice) %>%
#   summarise(
#     mean_RT = mean(last_RT, na.rm = TRUE),
#     sd_RT = sd(last_RT, na.rm = TRUE),
#     .groups = 'drop'
#   )
# 
# # Tracer un diagramme en barres des moyennes avec barres d'erreur
# ggplot(summary_data, aes(x = as.factor(last_binary_choice), y = mean_RT, fill = Session)) +
#   geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
#   geom_errorbar(aes(ymin = mean_RT - sd_RT, ymax = mean_RT + sd_RT),
#                 width = 0.2, position = position_dodge(width = 0.9)) +
#   labs(title = "Temps de Réaction Moyen par Choix Binaire",
#        x = "Choix Binaire", y = "Temps de Réaction Moyen") +
#   # scale_fill_brewer(palette = "Set1") +
#   theme_minimal()


############# âge des participants ###########
# Importer les bibliothèques nécessaires
library(ggplot2)

# Données d'âge des participants
ages <- c(23, 23, 24, 50, 19, 53, 23, 18, 23, 25, 22)

# Calculer la moyenne et l'écart-type
mean_age <- mean(ages)
sd_age <- sd(ages)

# Créer un dataframe
data <- data.frame(Age = ages)

# Tracer un diagramme en violon
ggplot(data, aes(x = "", y = Age)) +
  geom_violin(fill = "skyblue", color = "blue") +
  geom_point(aes(y = mean_age), color = "red", size = 3, shape = 18) +
  geom_errorbar(aes(ymin = mean_age - sd_age, ymax = mean_age + sd_age), 
                color = "red", width = 0.2, size = 1) +
  geom_hline(yintercept = mean_age, linetype = "dashed", color = "red") +
  labs(title = "Distribution des âges des participants",
       y = "Âge",
       x = "") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


#################### STPQ des participants ####################

library(dplyr)
library(tidyr)
library(ggplot2)

# Données des sujets
data <- data.frame(
  Subject = paste0("S", 1:12),
  Total_T_Score = c(58, 63, 44,51,48, 54, 55, 47, 51, 57, 54, 35),
  F1_T_Score = c(60, 56, 38,64,49, 50, 52, 48, 50, 55, 46, 33),
  F2_T_Score = c(62, 67,49, 38, 55,57, 58, 60, 53, 58, 62, 35),
  F3_T_Score = c(43, 58, 53,50,38, 51, 51, 32, 48, 53, 53, 50)
)

# Calcul de la moyenne et de l'écart-type
mean_score <- 50
sd_score <- 10

# Conversion des données du format large au format long
data_long <- data %>%
  pivot_longer(cols = starts_with("F"), names_to = "Factor", values_to = "T_Score")

# Groupement
data_long <- data_long %>%
  mutate(Group = ifelse(T_Score >= mean_score, "Au-dessus de la moyenne", "En-dessous de la moyenne"))

# Ajouter Total_T_Score aux données au format long
total_scores <- data %>%
  select(Subject, Total_T_Score) %>%
  pivot_longer(cols = Total_T_Score, names_to = "Factor", values_to = "T_Score") %>%
  mutate(Group = ifelse(T_Score >= mean_score, "Au-dessus de la moyenne", "En-dessous de la moyenne"))

# Fusion des données
data_long <- bind_rows(data_long, total_scores)

# Tracer un diagramme en barres
ggplot(data_long, aes(x = Subject, y = T_Score, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~ Factor) +
  geom_hline(yintercept = mean_score, linetype = "dashed", color = "blue") +
  geom_errorbar(aes(ymin = mean_score - sd_score, ymax = mean_score + sd_score), 
                color = "red", width = 0.2, size = 1) +
  labs(title = "Scores T des participants au test STPQ",
       y = "Score T",
       x = "Participant") +
  theme_minimal() +
  scale_fill_manual(values = c("En-dessous de la moyenne" = "skyblue", "Au-dessus de la moyenne" = "orange"))

################# Lecture des données HHMT ####################

# Créer un dataframe
data <- data.frame(
  Participant = c("P08", "P09", "P10", "P04", "P05", "P11", 
                  "P07", "P13", "P12", "P06", "P02", "P03", "P01"),
  TxRep = c(0.4878048780487805, 0.5365853658536586, 0.6341463414634146, 0.4634146341463415,
            0.7317073170731707, 0.4878048780487805, 0.4878048780487805, 0.8048780487804879,
            0.7804878048780488, 0.6829268292682927, 0.4878048780487805, 0.6585365853658537,
            0.5121951219512195)
)

# Tracer le diagramme de distribution
ggplot(data, aes(x = as.factor(1), y = TxRep)) +
  geom_violin(fill = "skyblue", color = "blue") +
  labs(title = "Distribution des taux de succès dans le test HHMT",
       y = "Taux de succès",
       x = "Participant") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

########################### Lecture des données HMMT ###########################

good_FilePath="/Users/wenjiehuang/Desktop/LyonLumiere2/M2/StageM2/analyse/comp/HHMT_good"
good_csv_files <- list.files(good_FilePath, pattern = "\\.csv$", full.names = TRUE)

# Initialiser un dataframe vide
good_data <- data.frame()

# Définir le taux d'échantillonnage original et le taux cible
original_sampling_rate <- 0.002
target_sampling_rate <- 0.05

# Calculer le nombre de lignes à conserver
rows_to_keep <- round(target_sampling_rate / original_sampling_rate)

# Parcourir chaque fichier csv
for (file in good_csv_files) {
  # Extraire le numéro du sujet et le scénario de l'expérience à partir du nom de fichier
  file_name <- basename(file)
  subject_number <- strsplit(file_name, "_")[[1]][1]
  scenario <- strsplit(file_name, "_")[[1]][2]
  
  # Lire le fichier csv
  csv_data <- read.csv(file, header = TRUE, sep = ";", dec = ",")
  
  # Ajouter les colonnes subject_number et session
  csv_data$Subject <- subject_number
  
  # Échantillonner les données
  sampled_data <- csv_data[seq(1, nrow(csv_data), by = rows_to_keep), ]
  
  # Assurer que TypeCond et NbTrial sont du bon type
  sampled_data$TypeCond <- as.character(sampled_data$TypeCond)
  sampled_data$NbTrial <- as.numeric(sampled_data$NbTrial)
  
  # Ajouter TrialOrder en fonction du scénario et du numéro du sujet
  sampled_data <- sampled_data %>%
    group_by(Subject) %>%
    arrange(Time) %>%
    mutate(
      TrialOrder = 1 + cumsum(if_else(TypeCond != lag(TypeCond, default = first(TypeCond)) |
                                        NbTrial != lag(NbTrial, default = first(NbTrial)), 1, 0))
    )%>% 
    ungroup()
  
  # Ajouter les données au dataframe total
  good_data <- bind_rows(good_data, sampled_data)
}

# Remplacer toutes les valeurs 2 dans BinaryChoice par 0.5
good_data <- good_data %>%
  mutate(BinaryChoice = ifelse(BinaryChoice == 2, 0.5, BinaryChoice))

# Conversion initiale
good_data_modified <- good_data %>%
  mutate(
    BinaryChoice = 2 * as.numeric(BinaryChoice) - 1,
    JoystickYAxis = 2 * JoystickYAxis - 1
  )

# Modifier les données
good_data_modified <- good_data_modified %>%
  group_by(Subject, TrialOrder) %>%
  mutate(
    last_binary_choice = last(BinaryChoice),
    last_RT=last(RT)
  ) %>%
  filter(last_binary_choice != 0)

# Supprimer les lignes commençant par "SameSpeed" et marquer "SpeedUp" et "SlowDown"
good_data_modified <- good_data_modified %>%
  filter(!startsWith(TypeCond, "SameSpeed")) %>%
  mutate(TypeCond = case_when(
    startsWith(TypeCond, "SpeedUp") ~ "SpeedUp",
    startsWith(TypeCond, "SlowDown") ~ "SlowDown",
    TRUE ~ TypeCond
  ))

# Calcul du taux de succès
goodHH_HMMT_success_rate <- good_data_modified %>%
  mutate(correct = case_when(
    TypeCond == "SlowDown" & last_binary_choice == -1 ~ 1,
    TypeCond == "SpeedUp" & last_binary_choice == 1 ~ 1,
    TRUE ~ 0
  )) %>%
  group_by(Subject) %>%
  summarise(HMMT_success_rate = mean(correct))

# Afficher les résultats
print(goodHH_HMMT_success_rate)

########################### Lecture des données HMMT pour badHH ###########################

bad_FilePath="/Users/wenjiehuang/Desktop/LyonLumiere2/M2/StageM2/analyse/comp/HHMT_bad"
bad_csv_files <- list.files(bad_FilePath, pattern = "\\.csv$", full.names = TRUE)

# Initialiser un dataframe vide
bad_data <- data.frame()

# Définir le taux d'échantillonnage original et le taux cible
original_sampling_rate <- 0.002
target_sampling_rate <- 0.05

# Calculer le nombre de lignes à conserver
rows_to_keep <- round(target_sampling_rate / original_sampling_rate)

# Parcourir chaque fichier csv
for (file in bad_csv_files) {
  # Extraire le numéro du sujet et le scénario de l'expérience à partir du nom de fichier
  file_name <- basename(file)
  subject_number <- strsplit(file_name, "_")[[1]][1]
  
  # Lire le fichier csv
  csv_data <- read.csv(file, header = TRUE, sep = ";", dec = ",")
  
  # Ajouter les colonnes subject_number et session
  csv_data$Subject <- subject_number
  
  # Échantillonner les données
  sampled_data <- csv_data[seq(1, nrow(csv_data), by = rows_to_keep), ]
  
  # Assurer que TypeCond et Nb
  
  Trial sont du bon type
  sampled_data$TypeCond <- as.character(sampled_data$TypeCond)
  sampled_data$NbTrial <- as.numeric(sampled_data$NbTrial)
  
  # Ajouter TrialOrder en fonction du scénario et du numéro du sujet
  sampled_data <- sampled_data %>%
    group_by(Subject) %>%
    arrange(Time) %>%
    mutate(
      TrialOrder = 1 + cumsum(if_else(TypeCond != lag(TypeCond, default = first(TypeCond)) |
                                        NbTrial != lag(NbTrial, default = first(NbTrial)), 1, 0))
    )%>% 
    ungroup()
  
  # Ajouter les données au dataframe total
  bad_data <- bind_rows(bad_data, sampled_data)
}

# Remplacer toutes les valeurs 2 dans BinaryChoice par 0.5
bad_data <- bad_data %>%
  mutate(BinaryChoice = ifelse(BinaryChoice == 2, 0.5, BinaryChoice))

# Conversion initiale
bad_data_modified <- bad_data %>%
  mutate(
    BinaryChoice = 2 * as.numeric(BinaryChoice) - 1,
    JoystickYAxis = 2 * JoystickYAxis - 1
  )

# Modifier les données
bad_data_modified <- bad_data_modified %>%
  group_by(Subject, TrialOrder) %>%
  mutate(
    last_binary_choice = last(BinaryChoice),
    last_RT=last(RT)
  ) %>%
  filter(last_binary_choice != 0)

# Supprimer les lignes commençant par "SameSpeed" et marquer "SpeedUp" et "SlowDown"
bad_data_modified <- bad_data_modified %>%
  filter(!startsWith(TypeCond, "SameSpeed")) %>%
  mutate(TypeCond = case_when(
    startsWith(TypeCond, "SpeedUp") ~ "SpeedUp",
    startsWith(TypeCond, "SlowDown") ~ "SlowDown",
    TRUE ~ TypeCond
  ))

# Calcul du taux de succès
badHH_HMMT_success_rate <- bad_data_modified %>%
  mutate(correct = case_when(
    TypeCond == "SlowDown" & last_binary_choice == -1 ~ 1,
    TypeCond == "SpeedUp" & last_binary_choice == 1 ~ 1,
    TRUE ~ 0
  )) %>%
  group_by(Subject) %>%
  summarise(HMMT_success_rate = mean(correct))

# Afficher les résultats
print(badHH_HMMT_success_rate)

# Calcul du taux de succès global pour HH et HMMT
all_HMMT_success_rate <- bind_rows(
  goodHH_HMMT_success_rate %>% mutate(Group = "goodHH"),
  badHH_HMMT_success_rate %>% mutate(Group = "badHH")
)

# Afficher le taux de succès global
print(all_HMMT_success_rate)

# Tracer un diagramme de distribution
ggplot(all_HMMT_success_rate, aes(x = Group, y = HMMT_success_rate)) +
  geom_violin(fill = "skyblue", color = "blue") +
  labs(title = "Distribution des taux de succès dans le test HMMT par groupe",
       y = "Taux de succès",
       x = "Groupe") +
  theme_minimal()
