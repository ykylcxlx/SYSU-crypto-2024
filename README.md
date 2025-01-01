# SYSU-crypto-2024

实验网站为[现代密码学实验](https://crypto-lab.akarin.dev/#/)，可见各次实验要求

因为实验网站只提交代码允许提交一个文件，所有每次实验的内容都写在一个.cpp文件里面，显得很丑陋

所有的代码都通过了oj测评，AES和SHA-256分别有普通版本和NI指令集加速实现。

大数运算部分参考了[大数运算实现 | Smallorange&#39;s Blog](https://smallorange666.github.io/2024/11/29/%E5%A4%A7%E6%95%B0%E8%BF%90%E7%AE%97%E5%AE%9E%E7%8E%B0/)，使用了蒙哥马利模幂加速

持续更新中

TODO:

sha-256没有实现init函数，每次调用final之后要自行重置一下bitlen等变量。
