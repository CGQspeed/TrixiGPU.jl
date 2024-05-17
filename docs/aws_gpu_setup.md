# Tutorial - Setting Up Cloud GPUs for Running CUDA.jl

Welcome to this tutorial! If you're looking to use GPUs in the cloud for running Julia programs with CUDA.jl, you're in the right place. This tutorial aims to help you set up cloud GPU on Amazon Web Service (AWS) to run CUDA.jl.

## Prerequisite

Before the tutorial begins, you should launch P type EC2 instance(s) on your own in AWS. Here are some links for you to quickly get to know key concepts and suceessfully launch your instance(s).

- What is Amazon EC2? https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/concepts.html
- Which type of EC2 should I use?
  https://aws.amazon.com/ec2/instance-types/
- How to launch Amazon EC2? https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/EC2_GetStarted.html
- What is spot instance? (optional) https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/using-spot-instances.html
- How to work with spot instance? (optional) https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/spot-requests.html

The last two resources are optional if the cost isn't a concern for you, and you can choose to stick with on-demand instance(s). However, please note that in both cases you will need to manually request a quota for P type instances before launching.

## Connect to Instance(s)

Up to this step, we assume that you have already launched your EC2 instance(s). Now, you can connect to your instance(s) via your local termial or any IDE with an SSH extension.

Move your key pair file (.pem) to the `.ssh` folder:

```Bash
$ mv ~/Downloads/<your-key-pair-name>.pem ~/.ssh/<your-key-pair-name>.pem
```

Change the permission of the key pair file:

```Bash
$ chmod 600 ~/.ssh/<your-key-pair-name>.pem
```

Copy the public IPv4 DNS of your instance(s) from AWS and then connect to your instance(s) via SSH:

```Bash
$ ssh -i ~/.ssh/<your-key-pair-name>.pem ubuntu@<your-ec2-public-ipv4-dns>
```

When prompted with the question `Are you sure you want to continue connecting (yes/no)?`, enter `yes`. At this point, you should be able to successfully connect to your instance(s).

## Configure Julia Package

In this step, we assume that you have already connected to your instance(s) and the commands are being executed in your instance(s)'s terminal.

You can find the version of the Julia package you want from https://julialang.org/downloads/. For illustration, we are using Julia version 1.10.0 for Linux on x86 64-bit.

Copy the address of the Julia package and download it to your instance(s):

```Bash
$ wget https://julialang-s3.julialang.org/bin/linux/x64/1.10/julia-1.10.0-linux-x86_64.tar.gz
```

Extract your downloaded package:

```Bash
$ tar -xvf julia-1.10.0-linux-x86_64.tar.gz
```

Create a symbolic link to the Julia:

```Bash
$ sudo ln -s ~/julia-1.10.0/bin/julia /usr/local/bin/julia
```

Remove the original Julia package (optional):

```Bash
$ rm julia-1.10.0-linux-x86_64.tar.gz
```

Then you should be able to run `julia` command successfully in your instance(s)'s root directory.

## Configure CUDA Toolkit

Like the former step, we assume that you have connected to your instance(s) and the commands are being executed in your instance(s)'s terminal.

Then you can easily download and configure the desired version of the CUDA Toolkit based on your selected platform through a set of commands available at https://developer.nvidia.com/cuda-downloads. After that, you may verify the version of the CUDA Toolkit you installed by running the `nvcc --version` command.

## Add CUDA Toolkit to Julia

In this step, we are going to add CUDA package to Julia. Again, this is based on the fact that you have already connected to your insrance(s).

Enter into Julia REPL:

```Julia
$ julia
```

Enter Jualia package mode and add CUDA to Julia:

```Julia
pkg> add CUDA
```

Then you can simply test CUDA in Julia by creating a `test.jl` file like below:

```Julia
# This is a test.jl file
using CUDA
CUDA.versioninfo()
```

If the test result displays the version of CUDA without any errors, you have successfully added CUDA to Julia.

## Enable SSH with Git Repository (Optional)

This step is optional but recommended for those who are going to operate their git repository on EC2 instance(s).

Create an SSH key in your local terminal:

```Bash
$ ssh-keygen -t rsa -b 4096 -C <your-email-address>
```

Save the key to the default directory `~/.ssh` and enter a passphrase of your choice.

Change the permission of your SSH key:

```Bash
$ chmod 600 ~/.ssh/id_rsa
```

Add your SSH private key to the ssh-agent:

```Bash
$ ssh-add -k ~/.ssh/id_rsa
```

Copy the public key and use it to create a new SSH key in your GitHub.

Copy the public key into your EC2 instance(s):

```Bash
$ cat ~/.ssh/id_rsa.pub | ssh -i ~/.ssh/<your-key-pair-name>.pem ubuntu@<your-ec2-public-ipv4-dns> "cat >> .ssh/authorized_keys"
```

Copy the private key to your EC2 instance(s):

```Bash
$ scp -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no -i ~/.ssh/<your-key-pair-name>.pem ~/.ssh/id_rsa ubuntu@<your-ec2-public-ipv4-dns>:~/.ssh/
```

Then you can directly `git clone` and `git push` (or other operations) in your EC2 instance(s).

## In The End

After going through the previous steps, you can use CUDA.jl on your cloud GPU through AWS. Please remember to terminate your instance(s) when you no longer need them, as they will continue to incur charges.
