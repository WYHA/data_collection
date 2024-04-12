# 背景

天文学家通常对一些天文相关的python包有使用需求，比如根据条件进行相关数据检索、数据分析、表格分析与生成、读写文件、画图等。

**目标**：支持用户用语言表述的方式使用python包，增强在天文python库根据任务指令生成代码的能力。

**例子**

1. 输入用户任务指令：查询glimpse_s07星表中以给定txt文件中天体坐标为中心的2角分方形区域内的所有红外天体，txt地址在'/home/usr/object_coordinates.txt'，每一行包含一个天体坐标。

2. 根据用户的任务指令和依赖的python包生成python_code：
```python
def canonical_solution():
    from astroquery.ipac.irsa import Irsa
    import astropy.units as u
    Irsa.ROW_LIMIT = 10000
    try:
        with open('/home/usr/object_coordinates.txt', 'r') as file:
            object_coordinates_list = file.readlines()
    except:
        return 'File not found.'
    table = Irsa.query_region(object_coordinates_list, catalog="glimpse_s07", spatial='Box', width=2*u.arcmin).to_pandas()
    return table
```

3. 执行python_code返回table
![image.png](https://cdn.nlark.com/yuque/0/2024/png/29422557/1710900443535-d2185433-0352-4734-a368-7c34a97d16f0.png#averageHue=%23282828&clientId=u153c5ad7-f01b-4&from=paste&height=839&id=uf31b0723&originHeight=1258&originWidth=3155&originalType=binary&ratio=1.5&rotation=0&showTitle=false&size=278964&status=done&style=none&taskId=u01f11296-40e2-4df5-b411-db88997ea7c&title=&width=2103.3333333333335)

# 数据收集需求

我们希望收集一些任务需求制作上述任务的评测集，先从几个常用的天文python包入手：astropy，astroquery（如果还有别的常用python包欢迎补充）。需要收集的内容包括：**用户任务指令、依赖的python包、Python code、需求来源**。

* 用户任务指令：用自然语言的方式描述需要astropy，astroquery等完成的任务，组织成一个明确的文本（比如数据检索任务，文本应包含所有的检索条件），如有依赖的文件，可以上传至 input_data 下，并记录上传文件地址。例如：查询glimpse_s07星表中以给定txt文件中天体坐标为中心的2角分方形区域内的所有红外天体，txt地址在'/home/usr/object_coordinates.txt'，每一行包含一个天体坐标。

* 依赖的python包：完成任务需要的python包。例如：astroquery、astropy。

* python code：用来完成用户任务的python code。由于太长的python code很难用文本描述清楚，最好将python code的长度控制在100行以内。请将可以执行成功的python code填写在canonical_solution函数里，结果通过return返回。<br>"def canonical_solution():\n  ...(此处省略完成任务对应的python code)\n  return data"，直接将代码作为文本进行读取，将读取的文本内容复制粘贴至此列。

* 需求来源：添加数据的人员、来源链接等，请务必填写

| 用户任务指令 | 依赖的python包 | Python code | 需求来源 |
| --- | --- | --- | --- |
| 查询glimpse红外巡天星表中以给定txt文件中天体坐标为中心的2角分方形区域内的所有红外天体，txt地址在'/home/usr/object_coordinates.txt'，每一行包含一个天体坐标<br>**依赖文件上传至**：input_data/object_coordinates.txt | astropy、astroquery |"def canonical_solution():\n  from astroquery.ipac.irsa import Irsa\n  import astropy.units as u\n  Irsa.ROW_LIMIT = 10000\n  try:\n    with open('/home/usr/object_coordinates.txt', 'r') as file:\n      object_coordinates_list = file.readlines()\n  except:\n      return 'File not found.'\n  table = Irsa.query_region(object_coordinates_list, catalog='glimpse_s07', spatial='Box', width=2*u.arcmin).to_pandas()\n  return table\n"| 张天惟 |
| 查询glimpse红外巡天星表中以13:16:43.64 -62:58:31.39坐标为中心的2角分方形区域内的所有红外天体 | astroquery | "def canonical_solution():\n  from astroquery.ipac.irsa import Irsa\n  import astropy.units as u\n  Irsa.ROW_LIMIT = 10000\n  table = Irsa.query_region("13:16:43.64 -62:58:31.39", catalog="glimpse_s07", spatial='Box', width=2*u.arcmin).to_pandas()\n  return table\n" | 张天惟 |
| 将赤道坐标10.625，41.2转换为银道坐标 | astropy | "def canonical_solution():\n  from astropy import units as u\n  from astropy.coordinates import SkyCoord\n  c = SkyCoord(ra=10.625*u.degree, dec=41.2*u.degree, frame='icrs')\n  return c.galactic\n"| 张天惟<br>[https://docs.astropy.org/en/stable/coordinates/](https://docs.astropy.org/en/stable/coordinates/) |
| 将赤道坐标10.625，41.2转换为时分秒写法 | astropy | "def canonical_solution():\n  from astropy import units as u\n  from astropy.coordinates import SkyCoord\n  c = SkyCoord(ra=10.625*u.degree, dec=41.2*u.degree, frame='icrs')\n  c.to_string('hmsdms')\n  return c\n" |  张天惟|
| 计算ICRS坐标系下的坐标10，9和FK5坐标系下的坐标11，10之间的空间距离 | astropy | "def canonical_solution():\n  from astropy import units as u\n  from astropy.coordinates import SkyCoord\n  c1 = SkyCoord(ra=10*u.degree, dec=9*u.degree, frame='icrs')\n  c2 = SkyCoord(ra=11*u.degree, dec=10*u.degree, frame='fk5')\n  return c1.separation(c2)" |张天惟  |
| 将1970年的赤道坐标02h31m49.09s +89d15m50.8s修正为2000年的赤道坐标 | astropy | "def canonical_solution():\n  from astropy.time import Time\n  from astropy.coordinates import SkyCoord\n  from astropy.coordinates import FK5\n  fk5c = SkyCoord('02h31m49.09s', '+89d15m50.8s',               frame=FK5(equinox=Time('J1970')))\n  fk5_2000 = FK5(equinox=Time(2000, format='jyear'))\n  return fk5c.transform_to(fk5_2000)"| 张天惟 |
| 从'tutorials/FITS-images/HorseHead.fits'下载fits文件，做出灰度图和colorbar并保存为horsehead.png | astropy, matplotlib | "def canonical_solution():\n  import matplotlib.pyplot as plt\n  from astropy.visualization \nmport astropy_mpl_style\n  from astropy.io import fits\n  from astropy.utils.data import \net_pkg_data_filename\n  plt.style.use(astropy_mpl_style)\n  image_file = get_pkg_data_filename('tutorials/FITS-images/HorseHead.fits')\n  image_data = fits.getdata(image_file, ext=0)\n  plt.figure()\n  plt.imshow(image_data, cmap='gray')\n  plt.colorbar()\n  plt.savefig('horsehead.png')\n  return" |张天惟<br> [https://docs.astropy.org/en/stable/generated/examples/io/plot_fits-image.html](https://docs.astropy.org/en/stable/generated/examples/io/plot_fits-image.html) |
| 将地址为'HorseHead.fits'的fits文件做出灰度图，横轴纵轴为真实的天体坐标，加网格线，并保存为horsehead.png | astropy, matplotlib | "def canonical_solution():\n  from astropy.wcs import WCS\n  from astropy.io import fits\n  import matplotlib.pyplot as plt\n  hdu = fits.open('HorseHead.fits')[0]\n  wcs = WCS(hdu.header)\n  plt.subplot(projection=wcs) \n  plt.imshow(hdu.data, origin='lower', cmap='gray') \n  plt.grid(color='white', ls='solid')\n  plt.xlabel('RA (J2000)')\n  plt.ylabel('DEC (J2000)')\n  plt.savefig('horsehead.png')\n  return" | 张天惟 |
| 将地址为'HorseHead.fits'的fits文件做出灰度图，在最强发射80%的范围内按对数画出5条等高线，横轴纵轴为真实的天体坐标，并保存为horsehead.png | astropy, matplotlib, numpy | "def canonical_solution():\n  from astropy.wcs import WCS\n  from astropy.io import fits\n  import numpy as np\n  import matplotlib.pyplot as plt\n  hdu = fits.open('HorseHead.fits')[0]\n  wcs = WCS(hdu.header)\n  plt.subplot(projection=wcs) \n  plt.imshow(hdu.data, origin='lower', cmap='gray') \n  max_v = np.nanmax(hdu.data)\n  min_v = np.nanmin(hdu.data)\n  plt.contour(hdu.data,levels=np.logspace(min_v,max_v*0.8,5))\n  plt.xlabel('RA (J2000)')\n  plt.ylabel('DEC (J2000)')\n  plt.savefig('horsehead.png')\n  return" |张天惟<br> [https://astropy-astrofrog.readthedocs.io/en/latest/visualization/wcsaxes/images_contours.html](https://astropy-astrofrog.readthedocs.io/en/latest/visualization/wcsaxes/images_contours.html) |

可以直接在上表中新增行填写
