import setuptools  # 导入setuptools打包工具

#with open("README.md", "r", encoding="utf-8") as fh:
#    long_description = fh.read()

setuptools.setup(
    name="Ecoli_CRISPRtyping",
    version="1.0.0",  # 包版本号，便于维护版本,保证每次发布都是版本都是唯一的
    author="",  # 作者，可以写自己的姓名
    author_email="",  # 作者联系方式，可写自己的邮箱地址
    description="CRISPR typing of E.coli",  # 包的简述
    long_description='long_description',  # 包的详细介绍，一般在README.md文件内
    long_description_content_type="text/markdown",
    url="http",  # 自己项目地址，比如github的项目地址
    packages=setuptools.find_namespace_packages(),
    install_requires=['pandas','numpy','docopt','openpyxl'],
    entry_points={
          'console_scripts': [
              'Ecoli_CRISPRtyping= Ecoli_CRISPRtyping.core:main'
          ],

    },
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
