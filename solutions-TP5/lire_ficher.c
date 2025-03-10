#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#define MAX_LINE_LENGTH 1000
#define MAX_NUMBERS 1000

// Fonction pour vérifier si une ligne contient uniquement des commentaires ou est vide
int is_valid_line(char *line) {
    // Ignore les lignes vides
    if (line[0] == '\0' || line[0] == '\n') {
        return 0;
    }

    // Ignore les lignes de commentaires (on suppose que le commentaire commence par '#')
    if (line[0] == '#') {
        return 0;
    }

    return 1;
}

// Fonction pour extraire les nombres d'une ligne
int extract_numbers(char *line, int numbers[], int *count) {
    char *token;
    int num;
    int local_count = 0;

    // Séparer la ligne par les espaces
    token = strtok(line, " \t\n");

    // Parcourir chaque token et vérifier si c'est un nombre
    while (token != NULL) {
        // Vérifier si le token est un nombre entier
        if (sscanf(token, "%d", &num) == 1) {
            numbers[local_count++] = num;
        }
        token = strtok(NULL, " \t\n");
    }

    *count = local_count; // Retourner le nombre de nombres trouvés
    return local_count;
}

// Fonction pour convertir un nombre en binaire et l'enregistrer dans un tableau
void convert_to_binary(int number, char binary[32]) {
    for (int i = 31; i >= 0; i--) {
        binary[31 - i] = (number & (1 << i)) ? '1' : '0';
    }
    binary[32] = '\0'; // Ajouter le caractère de fin de chaîne
}

int main() {
    FILE *file = fopen("table5.txt", "r");  // Ouvrir le fichier
    if (file == NULL) {
        perror("Erreur d'ouverture du fichier");
        return 1;
    }

    char line[MAX_LINE_LENGTH];
    int numbers[MAX_NUMBERS]; // Tableau pour stocker les nombres
    int count = 0; // Compteur de nombres lus
    char binary[32]; // Tableau pour stocker la représentation binaire

    // Lire le fichier ligne par ligne
    while (fgets(line, MAX_LINE_LENGTH, file)) {
        // Vérifier si la ligne est valide (non vide et non un commentaire)
        if (is_valid_line(line)) {
            // Extraire les nombres de la ligne
            int local_count = 0;
            int numbers_in_line = extract_numbers(line, numbers, &local_count);
            if (numbers_in_line > 0) {
                // Pour chaque nombre extrait, le convertir en binaire et l'enregistrer
                for (int i = 0; i < numbers_in_line; i++) {
                    convert_to_binary(numbers[i], binary);
                    printf("Nombre: %d -> Binaire: %s\n", numbers[i], binary);
                }
                count += numbers_in_line; // Mettre à jour le compteur total
            }
        }
    }

    fclose(file); // Fermer le fichier

    printf("Nombre total de nombres lus: %d\n", count);
    return 0;
}
