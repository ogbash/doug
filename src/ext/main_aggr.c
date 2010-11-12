#include <stdio.h>

extern ext_doug_init_(int *init_type);

int main()
{
  int init_type=1;
  ext_doug_init_(&init_type); // parallel
  printf("Initialized\n");
  return 0;
}
