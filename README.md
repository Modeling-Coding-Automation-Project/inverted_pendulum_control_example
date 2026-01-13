# Furuta式（回転型）倒立振子の制御設計

Pythonで物理モデルを設計し、実行結果をmatplotlibで可視化できるようにしました。

また、制御モデルをPythonベースで設計し、CopilotでC++にコード生成できます。

生成したC++コードの振る舞いがPythonと一致しているかどうかを、SIL検証できます。

詳細については以下の記事をご参照ください。

https://note.com/claude_a/n/nef52249cd84a

## 環境構築手順

ベースの環境構築について、以下のリンク先ページに記載しています。

[環境構築手順](https://github.com/Modeling-Coding-Automation-Project/MCAP_repo_manager/blob/main/documents/environment.md)

上記リンク先を参考に、WSL Ubuntu 24.04 をインストールし、Dockerをインストールしてください。

その後、本リポジトリを以下のコマンドでクローンします。

``` bash
git clone https://github.com/Modeling-Coding-Automation-Project/inverted_pendulum_control_example.git
```

その後、サブモジュールの更新も必要なので、以下のコマンドで行います。

``` bash
cd ./inverted_pendulum_control_example/
git submodule update --progress --init
```

その後、inverted_pendulum_control_exampleディレクトリをVisual Studio Codeで開きます。

Visual Studio Codeの「Reopen in Container」を実行し、Docker環境を起動します。

## サポート
新規にissueを作成して、詳細をお知らせください。

## 貢献
コミュニティからのプルリクエストを歓迎します。もし大幅な変更を考えているのであれば、提案する修正についての議論を始めるために、issueを開くことから始めてください。

また、プルリクエストを提出する際には、関連するテストが必要に応じて更新または追加されていることを確認してください。

## ライセンス

[MIT License](./LICENSE.txt)