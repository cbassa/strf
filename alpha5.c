#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

// Fixed mapping
const char mapping[] = "0123456789ABCDEFGHJKLMNPQRSTUVWXYZ";

// Function to convert alpha5 format to number
int alpha5_to_number(const char *s)
{
    // Split components
    char first_char = s[0];
    int digits = 0;
    sscanf(s + 1, "%4d", &digits);  // Get the last four digits
    
    // Decode first character
    int first_digit = 0;
    if (first_char == ' ') {
        first_digit = 0;
    } else if (isdigit(first_char)) {
        first_digit = first_char - '0';
    } else {
        for (int i = 0; mapping[i] != '\0'; i++) {
            if (mapping[i] == first_char) {
                first_digit = i;
                break;
            }
        }
    }

    // Create and return number
    return first_digit * 10000 + digits;
}

// Function to convert number to alpha5 format
void number_to_alpha5(int number, char *result)
{
    // Split components
    int first_digit = number / 10000;
    int digits = number % 10000;

    // Get char
    char first_char = mapping[first_digit];

    // Create string
    sprintf(result, "%c%04d", first_char, digits);

    return;
}

// Function to strip leading spaces
void strip_leading_spaces(const char *s, char *result)
{
    while (*s == ' ') {
        s++;  // Skip leading spaces
    }
    strcpy(result, s);  // Copy the remainder of the string

    return;
}

// Function to zero pad a string
void zero_pad(const char *s, char *result)
{
    char stripped[6];  // Temporarily store the string without leading spaces
    strip_leading_spaces(s, stripped);

    if (isdigit(stripped[0])) {
        int value = atoi(stripped);
        sprintf(result, "%05d", value);
    } else {
        strcpy(result, stripped);  // No zero-padding for non-digit strings
    }
}

int oldmain() {
    const char *norad_strings[] = {"00001", "    1", " 1234", "05697", "12345", "A0000", "J2931", "W1928", "E8493", "P4018", "Z9999"};
    const int correct_numbers[] = {1, 1, 1234, 5697, 12345, 100000, 182931, 301928, 148493, 234018, 339999};
    
    for (int i = 0; i < 11; i++) {
        const char *norad_string = norad_strings[i];
        int correct_number = correct_numbers[i];

        int number = alpha5_to_number(norad_string);
        printf("|%s| -> |%6d|==|%6d| %d\n", norad_string, number, correct_number, number == correct_number);
    }

    printf("\n");

    for (int i = 0; i < 11; i++) {
        int correct_number = correct_numbers[i];
        char s[7];  // For storing the alpha5 string

        number_to_alpha5(correct_number, s);
        
        char zero_padded[6];  // For storing the zero-padded string
        zero_pad(norad_strings[i], zero_padded);

        printf("|%6d| -> |%s|==|%s| %d\n", correct_number, s, zero_padded, strcmp(s, zero_padded) == 0);
    }

    return 0;
}
