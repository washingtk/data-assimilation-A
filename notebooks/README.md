# notebooks

`課題名前.ipynb`等と分かりやすい名前でファイルを管理する。

>例
>> kadai1.ipynb
>> kadai1_runge_kutta.ipynb

`observed_data.pickle`は0.05(=６時間)ステップで、2年分積分したデータの後半1年分に、MTで生成した標準正規分布を足した「観測データ」。　　
```python  
with open("PATH_TO_observed_data.pickle", 'rb') as f:  
  dat = pickle.load(f)  
```
