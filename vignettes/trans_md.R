#>>>>> 函数说明 <<<<<<<<<<
# 将.Rmd文档转换为.md文档，从而快速在支持MD格式的博客发文
# 需要安装knitr包，命令为 install.packages("knitr")
# 参数说明：
# post_name: Rmd文件完整名称
# post_path: 生成md文件的保存路径，默认为当前工作路径

trans_md <- function(post_name=NULL,post_path=getwd()){
  # 注意函数中最好不要出现中文，否则可能报错
  if(is.null(post_name) | post_name==""){
    stop("A post name must be given!")
  }

  input_file   <- paste(post_path,post_name, sep="/")

  post_name <- gsub("\\.[Rr]md$", "", post_name)
  out_file     <- paste0(post_path, "/", post_name, ".md")

  knitr::knit(input = input_file, output = out_file)
  print("-------------DONE!--------------")
}
