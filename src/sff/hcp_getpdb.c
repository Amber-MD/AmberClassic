
/***********************************************************************
 *                          hcp_getpdb_prm()
 *
 * Appends to the tprmtop file.
 *
 * Calling Parameters: pdbfile (The input pdbfile)
 * 
 * return: None.
 ************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <ctype.h>
#include <time.h>


int main(int argc, char ** argv)
{
   FILE *fp_read,*fp_write;
   char str_buffer[85], tmp_str1[7], tmp_str2[5], tmp_str3[7], save_char;
   int i, strands=0, complex=0, residue=1, num_per_line;
   int n, res_start, num_blanks, skip_mult_blanks;

   if (argc != 3)
   {
      printf("ERROR: hcp_getpdb() requires 2 parameters - input pdb file name and output HCP prm file name\n");
      exit(1);
   }

   if ((fp_read = fopen(argv[1], "r")) == NULL)
   {
      printf("ERROR: hcp_getpdb() cannot open input pdb file");
      exit(1);
   }

   if((fp_write = fopen(argv[2], "w")) == NULL)
   {
      printf("ERROR: hcp_getpdb() cannot open output HCP prm file");
      exit(1);
   }

   fputs("%FLAG POINTERS FOR STRAND AND COMPLEX\n",fp_write);
   fputs("%FORMAT(10I8)\n",fp_write);

   /* count number of complexes and strands */
   while (!feof(fp_read)) 
   {
      fgets(str_buffer,72,fp_read);

      str_buffer[21]='\0';
      if(strcmp(str_buffer,"REMARK END-OF-COMPLEX") == 0)
      {
         complex++;
      }

      str_buffer[3]='\0';
      if(strcmp(str_buffer,"TER")==0)
      {
         strands++;
      }
   }

   /* in case of only one complex */
   if (!complex)
   {
      complex = 1;
   }


   /* write out number of complexes and strands */
   sprintf(tmp_str2,"%d",complex);
   for(i=0;i<8-strlen(tmp_str2);i++)
   {
      fputs(" ",fp_write);
   }
   fprintf(fp_write,"%d",complex);

   sprintf(tmp_str2,"%d",strands);
   for(i=0;i<8-strlen(tmp_str2);i++)
   {
      fputs(" ",fp_write);
   }
   fprintf(fp_write,"%d",strands);


   fputs("\n%FLAG STRAND_POINTER\n",fp_write);
   fputs("%FORMAT(10I8)\n",fp_write);
   fputs("       1",fp_write);
   tmp_str3[0]='\0';

   /* get starting residue number for each strand */
   rewind(fp_read);
   num_per_line = 0;
   n = 1;
   while ((!feof(fp_read)) && (n < strands)) 
   {
      fgets(str_buffer,72,fp_read);

      /* find start of chain or residue number */
      res_start = 0;
      num_blanks = 0;
      skip_mult_blanks = 0;
      for (i = 0; i < 30; i++)
      {
         if  ((str_buffer[i] == ' ') && (!(skip_mult_blanks)))
         {
            num_blanks++;
            skip_mult_blanks = 1;
         }
         else
         {
            if (str_buffer[i] != ' ')
            {
               skip_mult_blanks = 0;
               if (num_blanks == 4)
               {
                  res_start = i;
                  i = 30; /* done */
               }
            }
         }
      }
      if (res_start == 0)
      {
         printf("ERROR: start of residue/chain not found in %s \n", str_buffer);
         exit(1);
      }
             
      for (i = 0; i < 6; i++)
      {
         tmp_str1[i]=str_buffer[res_start+i];
      }
      tmp_str1[6]='\0';

      if(strcmp(tmp_str1,tmp_str3)!=0)
      {
         residue++;
      }
      strcpy(tmp_str3,tmp_str1);

      str_buffer[3]='\0';
      if(strcmp(str_buffer,"TER")==0)
      {
         n++;
         fgets(str_buffer,72,fp_read);
         fgets(str_buffer,72,fp_read);
         if(str_buffer[0]!='E')
         {			
            /* print 10 per line */
            num_per_line++;
            if (num_per_line >= 10)
            {
               fprintf(fp_write, "\n");
               num_per_line = 0;
            }

            sprintf(tmp_str2,"%d",residue);
            for(i=0;i<8-strlen(tmp_str2);i++)
            {
               fputs(" ",fp_write);
            }
            fprintf(fp_write,"%d",residue);
         }
      }
   }

   /* get starting strand for each complex */ 
   fputs("\n%FLAG COMPLEX_POINTER\n",fp_write);
   fputs("%FORMAT(10I8)\n",fp_write);
   fputs("       1",fp_write);
   rewind(fp_read);
   strands = 0;
   num_per_line = 0;
   n = 1;
   str_buffer[0] = '\0';
   while (!feof(fp_read) && (n < complex)) 
   {
      fgets(str_buffer,72,fp_read);
      save_char = str_buffer[3];
      str_buffer[3]='\0';
      if(strcmp(str_buffer,"TER")==0)
      {
      	strands++;
      }
      str_buffer[3] = save_char;
      str_buffer[21] = '\0';
      if (strcmp(str_buffer, "REMARK END-OF-COMPLEX") == 0)
      {
         n++;
         fgets(str_buffer,72,fp_read);
         str_buffer[3]='\0';
         if (strcmp(str_buffer,"END")!=0)
         {
            /* print 10 per line */
            num_per_line++;
            if (num_per_line >= 10)
            {
               fprintf(fp_write, "\n");
               num_per_line = 0;
            }

            sprintf(tmp_str2,"%d",strands+1);
      	    for(i=0;i<8-strlen(tmp_str2);i++)
            {
      	       fputs(" ",fp_write);
            }
   	       fprintf(fp_write,"%d",strands+1);
         }
      }
   }
   fprintf(fp_write, "\n");

   fclose(fp_read);
   fclose(fp_write);

   return(0);
}
