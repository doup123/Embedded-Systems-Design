#include <stdio.h>
#include <linux/kernel.h>
#include <sys/syscall.h>
#include <unistd.h>
 
#define __NR_hello 378 
 
long hello_syscall(void)
{
    return syscall(__NR_hello);
}
 
int main(int argc, char *argv[])
{
    long int a = hello_syscall();
    printf("System call returned %ld\n", a);
    return 0;
}
