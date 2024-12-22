#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int alpha5_to_number(char *s)
{
  int i;
  char c,spart[]="1234";
  char alpha5[]="ABCDEFGHJKLMNPQRSTUVWXYZ";
  
  printf("%c %s %d\n",alpha5[0],alpha5,strlen(alpha5));

  c = s[0];
  strncpy(spart,s+1,5);
  spart[5]='\0';
  for (i=0;i<strlen(alpha5);i++) 
    if (c==alpha5[i])
      break;

  printf(">> %d %c |%c| |%s|\n",i,c,alpha5[i],spart);
  return (i+10)*10000+atoi(spart);
}

void number_to_alpha5(int number,char *s)
{
  int i;
  char c,alpha5[]="ABCDEFGHJKLMNPQRSTUVWXYZ";

  i=number/10000;
  c=alpha5[i-10];
  sprintf(s,"%c%04d",c,number-i*10000);

  return;
}

int main(int argc, char *argv[])
{
  int number;
  char s[] = "A0000",st[]="     ";

  number = alpha5_to_number(s);

  printf("%d\n",number);

  number_to_alpha5(number,st);

  printf("%s\n",st);
  
  return 0;
}
