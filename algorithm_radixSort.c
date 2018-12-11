#include <stdio.h>
#define MAX 8
#define SHOWPASS

struct r {int key; float value;};

void print(struct r *a, int n)
{
  int i;
  for (i = 0; i < n; i++) {
    printf("{%d %03.0f}\t", a[i].key, a[i].value);
  }
}

void radixSortBinStruct(r *a, int n)
{
  int i;
  r q0[MAX], q1[MAX], b[MAX];
  int c0, c1; 
  int m = a[0].key; 
  int exp = 1;

  for (i = 0; i < n; i++)
  {
    if (a[i].key > m)
      m = a[i].key;
  }
 
  while (m / exp > 0)
  {
    c0=0;
    c1=0;
    for (i = 0; i < n; i++) {
      if((a[i].key / exp) % 2) {
	q0[c0++]=a[i];
      }else{
	q1[c1++]=a[i];
      }
    }
    if(c0+c1!=n) printf("ERROR c0+c1\n");

    for (i = 0; i < c0; i++)
      a[i]=q0[i];
    for (i = c0; i < n; i++)
      a[i]=q1[i-c0];

    exp *= 2; 
    #ifdef SHOWPASS
      printf("\nPASS   : ");
      print(a, n);
    #endif
  }
 
}

int main()
{
  r arr[MAX];
  int i, n;
 
  printf("Enter total elements (n < %d) : ", MAX);
  scanf("%d", &n);
 
  printf("Enter %d Elements : ", n);

  for (i = 0; i < n; i++)
    scanf("%d", &(arr[i].key));

  for (i = 0; i < n; i++)
    arr[i].value = (float) (arr[i].key * 10);
 
  printf("\nARRAY  : ");
  print(arr, n);
 
  radixSortBinStruct(arr, n);
 
  printf("\nSORTED : ");
  print(arr, n);
  printf("\n");
 
  return 0;
}
