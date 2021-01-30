# general_2d_NS_Chemical-N2-N-

How to use this code  
  
1 Put "general_2d_NS" anywhere you want  
  
2 Change the src_path in "general_2d_NS.jl"  
  
　　Specify the absolute path up to and including the src folder  
    
3 Change the number of cells etc. in "pre.jl" and run  
  
4 Specify the number of threads by running the following at the windows command prompt  
  
  However, due to the small number of grids, 2 or 1 is fine  
    
　　set JULIA_NUM_THREADS=1  
    
5 Run "general_2d_NS.jl"  
  
6 Run "post.jl"  
  
7 Post_result folder has the result, visualize them using paraview  
  
  
  <br>
  <br>
  <br>
    
コードの使い方  
  
1　general_2d_NSを好きなところに置く  
  
2　general_2d_NS.jl内のsrc_pathを変更する  
  
　　srcフォルダーの配下までを絶対パスで指定  
    
3　pre.jlのセル数等を変更し，実行  
  
4　windowsのコマンドプロンプトで下記を実行しスレッド数を指定  
  
  ただし，格子数が少ないため，2か1でよい  
    
　　set JULIA_NUM_THREADS=1  
    
5　general_2d_NS.jlを実行  
  
6　post.jlを実行  
  
7　post_resultフォルダーに結果があるので，paraviewを使って可視化  
  
  
This software is released under the MIT License, see LICENSE.
