# Aptamer
This repository is about the work that I'm trying to do to make aptamer selection easier and analysis better.

# Aptamer
迄今为止个人对于aptamer信息化、智能筛选做出的一些尝试和实践，共勉；

## **1.Aptamer正交组序列**

是对于一系列的相互正交的aptamer序列的主动生成的尝试；目的是基于一条长链，生成一系列与之完全正交的短aptamer序列来实现筛选设计、动态设计。

## **2.Aptamer-desgin**

是一系列尝试生成aptamer的pdb、mol格式文件以方便后续高通量、虚拟筛选的尝试。

## **3. 尝试使用linux异步调用实现Rosetta、pyrosetta、biopython的组合使用**   
在Linux系统中，分步调用不同的软件通常涉及到几个关键步骤：

  3.1 **安装软件**：首先确保你需要的软件已经安装在系统中。可以使用包管理器如`apt`、`yum`或`pacman`来安装软件。

3.2 **配置环境变量**：为了能够从任何位置调用软件，你需要将软件的可执行文件路径添加到环境变量`PATH`中。这可以通过编辑`~/.bashrc`或`~/.profile`文件来实现。例如：
   ```bash
   export PATH="/path/to/software:$PATH"
   ```
   然后，运行`source ~/.bashrc`来使改动生效¹。

3.3 **使用软件**：一旦软件路径被添加到`PATH`，你就可以直接通过软件名称来调用它，而不需要指定完整路径。

3.4 **进程间通信**：如果你需要让不同的软件之间进行通信，可以使用如命名管道（FIFO）、套接字（Socket）等多种进程间通信（IPC）机制²。

3.5 **编写脚本**：为了自动化这一过程，你可以编写一个脚本来顺序调用不同的软件，并处理它们之间的数据传递。

3.6 **权限管理**：确保你有足够的权限来执行这些软件。如果需要，可以使用`sudo`来获取超级用户权限。

3.7 **错误处理**：在脚本中添加错误处理机制，以确保在调用过程中能够捕获并处理可能出现的错误。


