#include <linux/kernel.h>

asmlinkage long sys_hello1(void)
{
    printk("Greeting from kernel and team no A17");
    return 0;
}
