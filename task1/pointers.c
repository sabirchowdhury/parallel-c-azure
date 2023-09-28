// Online C compiler to run C program online
#include <stdio.h>

typedef struct i {
    int *j;
    int *j2;
} i;

void swapJs (i *pop) {
    int *newj = pop->j;
    int *newj2 = pop->j2;
    
    int *temp = newj;
    newj = pop->j2;
    newj2 = pop->j;
    pop->j = newj;
    pop->j2 = newj2;

}

int main() {
    // Write C code here
    int x = 5;
    int y = 545;
    int *p = &x;

    i e = {.j = &x, .j2 = &y};
    i *l = &e;
    
    
    
    
    
    printf("Hello world j  %d\n",l->j);
    printf("Hello world j  %d\n",*l->j);
    printf("Hello world j2 %d\n",l->j2);
    printf("Hello world j2 %d\n",*l->j2);
    printf("Hello worle j  %d\n",e.j);
    printf("Hello worle j2 %d\n",e.j2);
    printf("Hello worle j  %d\n",*e.j);
    printf("Hello worle j2 %d\n",*e.j2);
    swapJs(l);
    printf("Hello world j  %d\n",l->j);
    printf("Hello world j  %d\n",*l->j);
    printf("Hello world j2 %d\n",l->j2);
    printf("Hello world j2 %d\n",*l->j2);
    printf("Hello worle j  %d\n",e.j);
    printf("Hello worle j2 %d\n",e.j2);
    printf("Hello worle j  %d\n",*e.j);
    printf("Hello worle j2 %d\n",*e.j2);

    return 0;
}