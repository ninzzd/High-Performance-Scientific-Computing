# include <math.h>
# include <unistd.h> 
# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include <time.h>

int main ( int argc, char *argv[] );
char ch_cap ( char ch );
int ch_eqi ( char ch1, char ch2 );
int ch_to_digit ( char ch );
void gmsh_data_read ( char *gmsh_filename, int node_dim, int node_num,
double node_x[], int element_order, int element_num, int element_node[] );
void gmsh_size_read ( char *gmsh_filename, int *node_num, int *node_dim,
int *element_num, int *element_order );
void i4mat_write ( char *output_filename, int m, int n, int table[] );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
void r8mat_write ( char *output_filename, int m, int n, double table[] );
int s_begin ( char *s1, char *s2 );
int s_len_trim ( char *s );
int s_to_i4 ( char *s, int *last, int *error );
double s_to_r8 ( char *s, int *lchar, int *error );
void timestamp ( );


int main ( int argc, char *argv[] )
{
    int *element_node;
    int element_num;
    int element_order;
    char fem_element_filename[512];
    char fem_node_filename[512];
    char gmsh_filename[512];
    int m;
    int node_num;
    double *node_x;
    char prefix[512];
    char cwd[512] = __FILE__;
    int null_index;
    // To find the parent directory of the source file
    for(int i = 511;i >= 0;i--){
        if(cwd[i] == '\0'){
            null_index = i;
            break;
        }
    }
    #ifdef _WIN32
    for(int j = null_index-1;j >= 0;j--){
        if(cwd[j] == '\\'){
            break;
        }
        else{
            cwd[j] = '\0';
        }
    }
    #elif __linux__
    for(int j = null_index-1;j >= 0;j--){
        if(cwd[j] == '/'){
            break;
        }
        else{
            cwd[j] = '\0';
        }
    }
    #endif
    printf("Source file directory: %s\n",cwd);

    timestamp ( );
    printf ( "\n" );
    if ( argc <= 1 )
    {
        printf ( "\n" );
        printf ( "  Please enter the filename prefix.\n" );
        scanf ( "%s", prefix );
    }
    else
    {
        strcpy ( prefix, argv[1] );
    }
    /*
    Create the filenames.
    */
    strcpy ( gmsh_filename, cwd );
    strcat ( gmsh_filename, prefix ); 
    strcat ( gmsh_filename, ".msh" );
    strcpy ( fem_node_filename, cwd );
    strcat ( fem_node_filename, prefix );
    strcat ( fem_node_filename, "_nodes.txt" );
    strcpy ( fem_element_filename, cwd );
    strcat ( fem_element_filename, prefix );
    strcat ( fem_element_filename, "_elements.txt" );
    printf("%s\n%s\n%s\n",gmsh_filename,fem_node_filename,fem_element_filename);
    /*
    Read GMSH sizes.
    */
    gmsh_size_read ( gmsh_filename, &node_num, &m, &element_num, &element_order );
    /*
    Report sizes.
    */
    printf ( "\n" );
    printf ( "  Size information from GMSH:\n" );
    printf ( "  Spatial dimension M = %d\n", m );
    printf ( "  Number of nodes NODE_NUM = %d\n", node_num );
    printf ( "  Number of elements ELEMENT_NUM = %d\n", element_num );
    printf ( "  Element order ELEMENT_ORDER = %d\n", element_order );
    FILE *f1;
    FILE *f2;int dilip;
    char nodeinfo[255];
    strcpy(nodeinfo,cwd);
    strcat(nodeinfo,"nodeinfo.txt");
    f1=fopen(nodeinfo,"w");
    for(dilip=0;dilip<2;dilip++)
    {
        if(dilip==0)
        {
            fprintf(f1, "%d\n",node_num );
        }
        if(dilip==1)
            fprintf(f1, "%d\n",m );
    }
    fclose(f1);
    char elementinfo[255];
    strcpy(elementinfo,cwd);
    strcat(elementinfo,"eleminfo.txt");
    f2=fopen(elementinfo,"w");
    for(dilip=0;dilip<2;dilip++)
    {
        if(dilip==0)
        {
            fprintf(f2, "%d\n",element_num );
        }
        if(dilip==1)
            fprintf(f2, "%d\n",element_order );
    }
    fclose(f2);


    /*
    Allocate memory.
    */
    node_x = ( double * ) malloc ( m * node_num * sizeof ( double ) );
    element_node = ( int * )
    malloc ( element_order * element_num * sizeof ( int ) );
    /*
    Read GMSH data.
    */
    gmsh_data_read ( gmsh_filename, m, node_num, node_x, element_order,
        element_num, element_node );
    /*
    Write FEM data.
    */
    r8mat_write ( fem_node_filename, m, node_num, node_x );

    i4mat_write ( fem_element_filename, element_order, element_num,
        element_node );
    /*
    Free memory.
    */
    free ( element_node );
    free ( node_x );
    /*
    Terminate.
    */
    printf ( "\n" );
    printf ( "GMSH_TO_FEM:\n" );
    printf ( "  Normal end of execution.\n" );
    printf ( "\n" );
    timestamp ( );

    return 0;
}
/******************************************************************************/

char ch_cap ( char ch )
{
  if ( 97 <= ch && ch <= 122 )
  {
    ch = ch - 32;
  }

  return ch;
}
/******************************************************************************/

int ch_eqi ( char ch1, char ch2 )

/******************************************************************************/

{
  int value;

  if ( 97 <= ch1 && ch1 <= 122 )
  {
    ch1 = ch1 - 32;
  }
  if ( 97 <= ch2 && ch2 <= 122 )
  {
    ch2 = ch2 - 32;
  }
  if ( ch1 == ch2 )
  {
    value = 1;
  }
  else
  {
    value = 0;
  }
  return value;
}
/******************************************************************************/

int ch_to_digit ( char ch )

/******************************************************************************/

{
  int digit;

  if ( '0' <= ch && ch <= '9' )
  {
    digit = ch - '0';
  }
  else if ( ch == ' ' )
  {
    digit = 0;
  }
  else
  {
    digit = -1;
  }

  return digit;
}
/******************************************************************************/

void gmsh_data_read ( char *gmsh_filename, int node_dim, int node_num,
  double node_x[], int element_order, int element_num, int element_node[] )

/******************************************************************************/

{
  char *error;
  int i;
  int ierror;
  FILE *input;
  int j;
  int k;
  int length;
  int level;
  char text[255];
  char* text_pointer;
  double x;

  input = fopen ( gmsh_filename, "rt" );

  if ( ! input )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "GMSH_DATA_READ - Fatal error!\n" );
    fprintf ( stderr, "  Could not open input file \"%s\"\n", gmsh_filename );
    exit ( 1 );
  }

  level = 0;

  for ( ; ; )
  {
    text_pointer = text;
    error = fgets ( text_pointer, 255, input );

    if ( !error )
    {
      break;
    }

    if ( level == 0 )
    {
      if ( s_begin ( text_pointer, "$Nodes" ) )
      {
        level = 1;
        j = 0;
      }
    }
    else if ( level == 1 )
    {
      s_to_i4 ( text_pointer, &length, &ierror );
      level = 2;
    }
    else if ( level == 2 )
    {
      if ( s_begin ( text_pointer, "$EndNodes" ) )
      {
        break;
      }
      else
      {
        s_to_i4 ( text_pointer, &length, &ierror );
        text_pointer = text_pointer + length;
        for ( i = 0; i < node_dim; i++ )
        {
          x = s_to_r8 ( text_pointer, &length, &ierror );
          text_pointer = text_pointer + length;
          node_x[i+j*node_dim] = x;
        }
        j = j + 1;
      }
    }
  }
/*
  Now read element information.
*/
  level = 0;

  for ( ; ; )
  {
    text_pointer = text;
    error = fgets ( text_pointer, 255, input );

    if ( !error )
    {
      break;
    }

    if ( level == 0 )
    {
      if ( s_begin ( text_pointer, "$Elements" ) )
      {
        level = 1;
        j = 0;
      }
    }
    else if ( level == 1 )
    {
      s_to_i4 ( text_pointer, &length, &ierror );
      level = 2;
    }
    else if ( level == 2 )
    {
      if ( s_begin ( text_pointer, "$EndElements" ) )
      {
        break;
      }
      else
      {
        for ( k = 1; k <= 5; k++ )
        {
          s_to_i4 ( text_pointer, &length, &ierror );
          text_pointer = text_pointer + length;
        }
        for ( i = 0; i < element_order; i++ )
        {
          k = s_to_i4 ( text_pointer, &length, &ierror );
          text_pointer = text_pointer + length;
          element_node[i+j*element_order] = k;
        }
        j = j + 1;
      }
    }
  }

  fclose ( input );

  return;
}
/******************************************************************************/

void gmsh_size_read ( char *gmsh_filename, int *node_num, int *node_dim,
  int *element_num, int *element_order )

/******************************************************************************/

{
  char *error;
  int ierror;
  FILE *input;
  int k;
  int length;
  int level;
  const double r8_big = 1.0E+30;
  char text[255];
  char* text_pointer;
  double x;
  double x_max;
  double x_min;
  double y;
  double y_max;
  double y_min;
  double z;
  double z_max;
  double z_min;

  *node_num = 0;
  *node_dim = 0;

  x_max = - r8_big;
  x_min = + r8_big;
  y_max = - r8_big;
  y_min = + r8_big;
  z_max = - r8_big;
  z_min = + r8_big;

  input = fopen ( gmsh_filename, "rt" );

  if ( ! input )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "GMSH_SIZE_READ - Fatal error!\n" );
    fprintf ( stderr, "  Could not open input file \"%s\"\n", gmsh_filename );
    exit ( 1 );
  }

  level = 0;

  for ( ; ; )
  {
    text_pointer = text;
    error = fgets ( text_pointer, 255, input );

    if ( !error )
    {
      break;
    }

    if ( level == 0 )
    {
      if ( s_begin ( text_pointer, "$Nodes" ) )
      {
        level = 1;
      }
    }
    else if ( level == 1 )
    {
      *node_num = s_to_i4 ( text_pointer, &length, &ierror );
      level = 2;
    }
    else if ( level == 2 )
    {
      if ( s_begin ( text_pointer, "$EndNodes" ) )
      {
        break;
      }
      else
      {
        s_to_i4 ( text_pointer, &length, &ierror );
        text_pointer = text_pointer + length;
        x = s_to_r8 ( text_pointer, &length, &ierror );
        x_min = r8_min ( x_min, x );
        x_max = r8_max ( x_max, x );
        text_pointer = text_pointer + length;
        y = s_to_r8 ( text_pointer, &length, &ierror );
        y_min = r8_min ( y_min, y );
        y_max = r8_max ( y_max, y );
        text_pointer = text_pointer + length;
        z = s_to_r8 ( text_pointer, &length, &ierror);
        text_pointer = text_pointer + length;
        z_min = r8_min ( z_min, z );
        z_max = r8_max ( z_max, z );
      }
    }
  }
/*
  Make a very simple guess as to the dimensionality of the data.
*/
  *node_dim = 3;
  if ( z_max == z_min )
  {
    *node_dim = 2;
    if ( y_max == y_min )
    {
      *node_dim = 1;
    }
  }
/*
  Now read element information.
*/
  level = 0;

  for ( ; ; )
  {
    text_pointer = text;
    error = fgets ( text_pointer, 255, input );

    if ( !error )
    {
      break;
    }

    if ( level == 0 )
    {
      if ( s_begin ( text_pointer, "$Elements" ) )
      {
        level = 1;
      }
    }
    else if ( level == 1 )
    {
      *element_num = s_to_i4 ( text_pointer, &length, &ierror );
      level = 2;
    }
    else if ( level == 2 )
    {
      if ( s_begin ( text_pointer, "$EndElements" ) )
      {
        break;
      }
      else
      {
        k = 0;
        for ( ; ; )
        {
          s_to_i4 ( text_pointer, &length, &ierror );
          text_pointer = text_pointer + length;
          if ( ierror != 0 )
          {
            break;
          }
          k = k + 1;
        }
        *element_order = k - 5;
        break;
      }
    }
  }

  fclose ( input );

  return;
}
/******************************************************************************/

void i4mat_write ( char *output_filename, int m, int n, int table[] )

/******************************************************************************/

{
  int i;
  int j;
  FILE *output;
/*
  Open the file.
*/
  output = fopen ( output_filename, "wt" );

  if ( !output )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "I4MAT_WRITE - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the file \"%s\".\n", output_filename );
    exit ( 1 );
  }
/*
  Write the data.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      fprintf ( output, "  %d", table[i+j*m] );
    }
    fprintf ( output, "\n" );
  }
/*
  Close the file.
*/
  fclose ( output );

  return;
}
/******************************************************************************/

double r8_max ( double x, double y )

/******************************************************************************/

{
  double value;

  if ( y < x )
  {
    value = x;
  }
  else
  {
    value = y;
  }
  return value;
}
/******************************************************************************/

double r8_min ( double x, double y )

/******************************************************************************/

{
  double value;

  if ( y < x )
  {
    value = y;
  }
  else
  {
    value = x;
  }
  return value;
}
/******************************************************************************/

void r8mat_write ( char *output_filename, int m, int n, double table[] )

/******************************************************************************/

{
  int i;
  int j;
  FILE *output;
/*
  Open the file.
*/
  output = fopen ( output_filename, "wt" );

  if ( !output )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8MAT_WRITE - Fatal error!\n" );
    fprintf ( stderr,
      "  Could not open the output file \"%s\".\n", output_filename );
    exit ( 1 );
  }
/*
  Write the data.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      fprintf ( output, "  %g", table[i+j*m] );
    }
    fprintf ( output, "\n" );
  }
/*
  Close the file.
*/
  fclose ( output );

  return;
}
/******************************************************************************/

int s_begin ( char *s1, char *s2 )

/******************************************************************************/

{
  int i;
  int n1;
  int n2;

  n1 = strlen ( s1 );
  n2 = strlen ( s2 );

  if ( n1 < n2 )
  {
    return 0;
  }

  for ( i = 0; i < n2; i++ )
  {
    if ( ch_cap ( s1[i] ) != ch_cap ( s2[i] ) )
    {
      return 0;
    }
  }
  return 1;
}
/******************************************************************************/

int s_len_trim ( char *s )

/******************************************************************************/

{
  int n;
  char *t;

  n = strlen ( s );
  t = s + strlen ( s ) - 1;

  while ( 0 < n )
  {
    if ( *t != ' ' && *t != '\n' )
    {
      return n;
    }
    t--;
    n--;
  }

  return n;
}
/******************************************************************************/

int s_to_i4 ( char *s, int *last, int *error )

/******************************************************************************/

{
  char c;
  int i;
  int isgn;
  int istate;
  int ival;

  *error = 0;
  istate = 0;
  isgn = 1;
  i = 0;
  ival = 0;

  while ( *s )
  {
    c = s[i];
    i = i + 1;
/*
  Haven't read anything.
*/
    if ( istate == 0 )
    {
      if ( c == ' ' )
      {
      }
      else if ( c == '-' )
      {
        istate = 1;
        isgn = -1;
      }
      else if ( c == '+' )
      {
        istate = 1;
        isgn = + 1;
      }
      else if ( '0' <= c && c <= '9' )
      {
        istate = 2;
        ival = c - '0';
      }
      else
      {
        *error = 1;
        return ival;
      }
    }
/*
  Have read the sign, expecting digits.
*/
    else if ( istate == 1 )
    {
      if ( c == ' ' )
      {
      }
      else if ( '0' <= c && c <= '9' )
      {
        istate = 2;
        ival = c - '0';
      }
      else
      {
        *error = 1;
        return ival;
      }
    }
/*
  Have read at least one digit, expecting more.
*/
    else if ( istate == 2 )
    {
      if ( '0' <= c && c <= '9' )
      {
        ival = 10 * (ival) + c - '0';
      }
      else
      {
        ival = isgn * ival;
        *last = i - 1;
        return ival;
      }

    }
  }
/*
  If we read all the characters in the string, see if we're OK.
*/
  if ( istate == 2 )
  {
    ival = isgn * ival;
    *last = s_len_trim ( s );
  }
  else
  {
    *error = 1;
    *last = 0;
  }

  return ival;
}
/******************************************************************************/

double s_to_r8 ( char *s, int *lchar, int *error )

/******************************************************************************/

{
  char c;
  int ihave;
  int isgn;
  int iterm;
  int jbot;
  int jsgn;
  int jtop;
  int nchar;
  int ndig;
  double r;
  double rbot;
  double rexp;
  double rtop;
  char TAB = 9;

  nchar = s_len_trim ( s );
  *error = 0;
  r = 0.0;
  *lchar = -1;
  isgn = 1;
  rtop = 0.0;
  rbot = 1.0;
  jsgn = 1;
  jtop = 0;
  jbot = 1;
  ihave = 1;
  iterm = 0;

  for ( ; ; )
  {
    c = s[*lchar+1];
    *lchar = *lchar + 1;
/*
  Blank or TAB character.
*/
    if ( c == ' ' || c == TAB )
    {
      if ( ihave == 2 )
      {
      }
      else if ( ihave == 6 || ihave == 7 )
      {
        iterm = 1;
      }
      else if ( 1 < ihave )
      {
        ihave = 11;
      }
    }
/*
  Comma.
*/
    else if ( c == ',' || c == ';' )
    {
      if ( ihave != 1 )
      {
        iterm = 1;
        ihave = 12;
        *lchar = *lchar + 1;
      }
    }
/*
  Minus sign.
*/
    else if ( c == '-' )
    {
      if ( ihave == 1 )
      {
        ihave = 2;
        isgn = -1;
      }
      else if ( ihave == 6 )
      {
        ihave = 7;
        jsgn = -1;
      }
      else
      {
        iterm = 1;
      }
    }
/*
  Plus sign.
*/
    else if ( c == '+' )
    {
      if ( ihave == 1 )
      {
        ihave = 2;
      }
      else if ( ihave == 6 )
      {
        ihave = 7;
      }
      else
      {
        iterm = 1;
      }
    }
/*
  Decimal point.
*/
    else if ( c == '.' )
    {
      if ( ihave < 4 )
      {
        ihave = 4;
      }
      else if ( 6 <= ihave && ihave <= 8 )
      {
        ihave = 9;
      }
      else
      {
        iterm = 1;
      }
    }
/*
  Exponent marker.
*/
    else if ( ch_eqi ( c, 'E' ) || ch_eqi ( c, 'D' ) )
    {
      if ( ihave < 6 )
      {
        ihave = 6;
      }
      else
      {
        iterm = 1;
      }
    }
/*
  Digit.
*/
    else if ( ihave < 11 && '0' <= c && c <= '9' )
    {
      if ( ihave <= 2 )
      {
        ihave = 3;
      }
      else if ( ihave == 4 )
      {
        ihave = 5;
      }
      else if ( ihave == 6 || ihave == 7 )
      {
        ihave = 8;
      }
      else if ( ihave == 9 )
      {
        ihave = 10;
      }

      ndig = ch_to_digit ( c );

      if ( ihave == 3 )
      {
        rtop = 10.0 * rtop + ( double ) ndig;
      }
      else if ( ihave == 5 )
      {
        rtop = 10.0 * rtop + ( double ) ndig;
        rbot = 10.0 * rbot;
      }
      else if ( ihave == 8 )
      {
        jtop = 10 * jtop + ndig;
      }
      else if ( ihave == 10 )
      {
        jtop = 10 * jtop + ndig;
        jbot = 10 * jbot;
      }
    }
/*
  Anything else is regarded as a terminator.
*/
    else
    {
      iterm = 1;
    }
/*
  If we haven't seen a terminator, and we haven't examined the
  entire string, go get the next character.
*/
    if ( iterm == 1 || nchar <= *lchar + 1 )
    {
      break;
    }

  }
/*
  If we haven't seen a terminator, and we have examined the
  entire string, then we're done, and LCHAR is equal to NCHAR.
*/
  if ( iterm != 1 && (*lchar) + 1 == nchar )
  {
    *lchar = nchar;
  }
/*
  Number seems to have terminated.  Have we got a legal number?
  Not if we terminated in states 1, 2, 6 or 7!
*/
  if ( ihave == 1 || ihave == 2 || ihave == 6 || ihave == 7 )
  {
    *error = 1;
    return r;
  }
/*
  Number seems OK.  Form it.

  We have had some trouble with input of the form 1.0E-312.
  For now, let's assume anything less than 1.0E-20 is zero.
*/
  if ( jtop == 0 )
  {
    rexp = 1.0;
  }
  else
  {
    if ( jbot == 1 )
    {
      if ( jsgn * jtop < -20 )
      {
        rexp = 0.0;
      }
      else
      {
        rexp = pow ( ( double ) 10.0, ( double ) ( jsgn * jtop ) );
      }
    }
    else
    {
      if ( jsgn * jtop < -20 * jbot )
      {
        rexp = 0.0;
      }
      else
      {
        rexp = jsgn * jtop;
        rexp = rexp / jbot;
        rexp = pow ( ( double ) 10.0, ( double ) rexp );
      }
    }

  }

  r = isgn * rexp * rtop / rbot;

  return r;
}
/******************************************************************************/

void timestamp ( )

/******************************************************************************/

{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
