void cubeOfOdd(int n) {
    int i = 1;
    while (i < n) {
        printf("%d\n", pow(i, 3));
        i += 2;
    }
}


int checkPrime(int n) {
    int status = 1;
    for (int i = 1; i < n/2; i++) {
        if (n % i == 0) {
            status = 0;
            break;
        }
    }
    return status;
}

void introToCS330(int n) {
    if ((n % 7 == 0) && (n % 3 == 0)) {
        printf("UAB CS 330\n");
    }
    else if (n % 7 == 0) {
        printf("UAB\n");
    }
    else if (n % 3 == 0) {
        printf("CS\n");
    }
    else if (checkPrime(n) && n != 3 && n != 7) {
        printf("Go Blazers\n");
    }
    else {
        printf("%d\n", pow(n, 3));
    }

}


#include <string.h>
void printHello(int n) { // This one is quite difficult. It will take some time to think about it.
    char str[100]; 
    // WARNING: the size of the array should be assigned in a better way. 
    //  I would try to find a way to relate int n to the number of 'HELLO's needed.
    strcpy(str, "0");

    char strtmp[5];
    for (int i = 0; i <= n; i++) {
        // Need if-else to test if power of 2 and use sprintf or strcpy to put value into strtmp.
        strcat(str, strtmp);
    }
}


int paintGallons(float length, float width, float height) {
    // I am unsure if this will work correctly in all cases because of the rounding up.
    //  But I believe it may since casting to int always rounds down.
    return (int) ((2 * (length * height) + 2 * (width * height) + length * width) / 400) + 1;
}


void grader(float avg_exams, float avg_hw, int attendance) {
    int pass = 1;
    if (attendance < 0) {
        pass = 0;
    }
    else if (avg_exams <= 70 || avg_hw <= 70) {
        pass = 0;
    }
    else if (avg_exams <= 85 && avg_hw <= 85) {
        pass = 0;
    }

    if (pass = 1) {
        printf("PASS\n");
    }
    else {
        printf("FAIL\n");
    }
}
