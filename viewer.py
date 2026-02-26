"""
output.bin ビューア — 2D CFD シミュレーション結果の可視化ツール

使い方:
    python viewer.py                  # カレントディレクトリの output.bin を開く
    python viewer.py path/to/file.bin # 指定パスのファイルを開く

操作:
    - フレームスライダー: フレーム（時刻）の切替
    - ラジオボタン: 表示フィールド切替 (u / v / speed / p)
    - ベクトルON/OFF: 速度ベクトルの表示切替
    - 密度スライダー: ベクトルの間引き間隔を調整
    - スケールスライダー: ベクトルの矢印サイズを調整
"""

import sys
import struct
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons, CheckButtons
from matplotlib.colors import Normalize
import matplotlib

# 日本語フォント対応
matplotlib.rcParams['font.family'] = ['MS Gothic', 'Yu Gothic', 'Meiryo', 'DejaVu Sans']


def read_header(f):
    """ヘッダを読み込み (zMax, xMax, dx, dz) を返す"""
    raw = f.read(4 + 4 + 8 + 8)
    if len(raw) < 24:
        raise ValueError("ヘッダの読み込みに失敗しました。ファイルが壊れています。")
    zMax, xMax = struct.unpack('ii', raw[:8])
    dx, dz = struct.unpack('dd', raw[8:24])
    return zMax, xMax, dx, dz


def build_frame_index(f, header_size, zMax, xMax):
    """全フレームのファイルオフセットと時刻のリストを構築する"""
    frame_data_size = 8 + 3 * (zMax * xMax * 8)
    offsets = []
    times = []

    f.seek(0, 2)
    file_size = f.tell()
    f.seek(header_size)

    pos = header_size
    while pos + frame_data_size <= file_size:
        f.seek(pos)
        t_raw = f.read(8)
        if len(t_raw) < 8:
            break
        t = struct.unpack('d', t_raw)[0]
        offsets.append(pos)
        times.append(t)
        pos += frame_data_size

    return offsets, times


def read_frame(f, offset, zMax, xMax):
    """指定オフセットのフレームを読み込み (t, u, v, p) を返す"""
    n = zMax * xMax
    f.seek(offset)
    t = struct.unpack('d', f.read(8))[0]
    u = np.fromfile(f, dtype=np.float64, count=n).reshape(zMax, xMax)
    v = np.fromfile(f, dtype=np.float64, count=n).reshape(zMax, xMax)
    p = np.fromfile(f, dtype=np.float64, count=n).reshape(zMax, xMax)
    return t, u, v, p


def compute_field(u, v, p, field_name):
    """選択されたフィールドのデータと表示ラベルを返す"""
    if field_name == 'u':
        return u, 'u (x方向速度) [m/s]'
    elif field_name == 'v':
        return v, 'v (z方向速度) [m/s]'
    elif field_name == 'speed':
        return np.sqrt(u**2 + v**2), '速さ |v| [m/s]'
    elif field_name == 'p':
        return p, '圧力 p [Pa]'
    return u, 'u'


CMAPS = {
    'u': 'coolwarm',
    'v': 'coolwarm',
    'speed': 'magma',
    'p': 'viridis',
}


def main():
    # ファイルパスの決定
    if len(sys.argv) > 1:
        filepath = sys.argv[1]
    else:
        filepath = 'output.bin'

    try:
        f = open(filepath, 'rb')
    except FileNotFoundError:
        print(f"エラー: ファイル '{filepath}' が見つかりません。")
        print("使い方: python viewer.py [output.bin のパス]")
        sys.exit(1)

    # ヘッダ読み込み
    zMax, xMax, dx, dz = read_header(f)
    header_size = 24
    print(f"格子: {zMax} x {xMax}, dx={dx:.6f}, dz={dz:.6f}")

    # フレームインデックス構築
    offsets, times = build_frame_index(f, header_size, zMax, xMax)
    num_frames = len(offsets)
    if num_frames == 0:
        print("エラー: フレームデータがありません。")
        f.close()
        sys.exit(1)
    print(f"フレーム数: {num_frames}  (t = {times[0]:.3f} ~ {times[-1]:.3f} s)")

    # 初期フレーム読み込み
    current_field = 'speed'
    show_vectors = False
    t, u, v, p = read_frame(f, offsets[0], zMax, xMax)
    data, label = compute_field(u, v, p, current_field)

    # 座標軸の生成
    x_coords = np.linspace(0, xMax * dx, xMax)
    z_coords = np.linspace(0, zMax * dz, zMax)

    # --- matplotlib セットアップ ---
    fig, ax = plt.subplots(figsize=(10, 9))
    plt.subplots_adjust(left=0.08, bottom=0.25, right=0.80, top=0.93)

    # 座標範囲
    extent = [0, xMax * dx, zMax * dz, 0]

    im = ax.imshow(data, extent=extent, cmap=CMAPS[current_field],
                   aspect='equal', interpolation='bilinear')
    ax.set_xlabel('x [m]')
    ax.set_ylabel('z [m]')
    title = ax.set_title(f'{label}  (t = {t:.3f} s)')

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(label)

    # NaN/Inf オーバーレイ画像 (RGBA)
    overlay_data = np.zeros((zMax, xMax, 4), dtype=np.float32)
    overlay_im = ax.imshow(overlay_data, extent=extent, aspect='equal',
                           interpolation='nearest', zorder=2)

    # quiver 用（初期状態は空）
    quiver_obj = [None]

    # --- ウィジェット ---

    # フレームスライダー
    ax_slider = plt.axes([0.15, 0.14, 0.50, 0.03])
    slider = Slider(ax_slider, 'Frame', 0, max(num_frames - 1, 1),
                    valinit=0, valstep=1, valfmt='%d')

    # ベクトル密度スライダー (間引き間隔: 小さい=密、大きい=疎)
    ax_density = plt.axes([0.15, 0.09, 0.50, 0.03])
    density_slider = Slider(ax_density, '密度', 1, 100,
                            valinit=5, valstep=1, valfmt='%d')

    # ベクトルスケールスライダー
    ax_scale = plt.axes([0.15, 0.04, 0.50, 0.03])
    scale_slider = Slider(ax_scale, 'スケール', 0.1, 50.0,
                          valinit=5.0, valfmt='%.1f')

    # フィールド切替ラジオボタン
    ax_radio = plt.axes([0.83, 0.45, 0.15, 0.2])
    radio = RadioButtons(ax_radio, ('u', 'v', 'speed', 'p'), active=2)

    # ベクトル表示ON/OFF チェックボックス
    ax_check = plt.axes([0.83, 0.35, 0.15, 0.08])
    check = CheckButtons(ax_check, ['ベクトル'], [show_vectors])

    def draw_quiver(u_data, v_data, step, scale):
        """速度ベクトルを描画する"""
        # 古い quiver を削除
        if quiver_obj[0] is not None:
            quiver_obj[0].remove()
            quiver_obj[0] = None

        if not show_vectors:
            return

        step = max(int(step), 1)
        # 間引いた格子点
        z_idx = np.arange(0, zMax, step)
        x_idx = np.arange(0, xMax, step)
        X, Z = np.meshgrid(x_coords[x_idx], z_coords[z_idx])
        U = u_data[np.ix_(z_idx, x_idx)]
        V = v_data[np.ix_(z_idx, x_idx)]

        quiver_obj[0] = ax.quiver(X, Z, U, V,
                                  color='white', alpha=0.8,
                                  scale=scale, scale_units='xy',
                                  width=0.002, headwidth=3,
                                  headlength=4, headaxislength=3)

    def update(val=None):
        nonlocal current_field
        frame_idx = int(slider.val)
        t, u, v, p = read_frame(f, offsets[frame_idx], zMax, xMax)

        # 診断情報をコンソールに出力
        nan_u, nan_v, nan_p = np.isnan(u).sum(), np.isnan(v).sum(), np.isnan(p).sum()
        inf_u, inf_v, inf_p = np.isinf(u).sum(), np.isinf(v).sum(), np.isinf(p).sum()
        print(f"--- Frame {frame_idx} (t={t:.4f}s) ---")
        print(f"  u:  min={np.nanmin(u):.6e}, max={np.nanmax(u):.6e}, NaN={nan_u}, Inf={inf_u}")
        print(f"  v:  min={np.nanmin(v):.6e}, max={np.nanmax(v):.6e}, NaN={nan_v}, Inf={inf_v}")
        print(f"  p:  min={np.nanmin(p):.6e}, max={np.nanmax(p):.6e}, NaN={nan_p}, Inf={inf_p}")

        data, label = compute_field(u, v, p, current_field)

        # NaN/Inf マスクを作成
        nan_mask = np.isnan(data)
        inf_mask = np.isinf(data)
        total_nan = nan_mask.sum()
        total_inf = inf_mask.sum()

        # NaN/Inf を除いた有限値のみでカラーレンジを決定
        finite_data = data[np.isfinite(data)]
        if len(finite_data) > 0:
            display_data = np.where(np.isfinite(data), data, 0.0)
            vmin, vmax = finite_data.min(), finite_data.max()
        else:
            display_data = np.zeros_like(data)
            vmin, vmax = -1.0, 1.0
        if vmin == vmax:
            vmin -= 0.5
            vmax += 0.5

        im.set_data(display_data)
        im.set_cmap(CMAPS[current_field])
        im.set_norm(Normalize(vmin=vmin, vmax=vmax))

        # NaN/Inf オーバーレイを更新
        overlay = np.zeros((zMax, xMax, 4), dtype=np.float32)
        # NaN → 赤 (半透明)
        overlay[nan_mask] = [1.0, 0.0, 0.0, 0.7]
        # Inf → 緑 (半透明)
        overlay[inf_mask] = [0.0, 1.0, 0.0, 0.7]
        overlay_im.set_data(overlay)

        # タイトルに NaN/Inf 情報を追加
        title_text = f'{label}  (t = {t:.3f} s)'
        if total_nan > 0 or total_inf > 0:
            title_text += f'  ⚠ NaN:{total_nan} Inf:{total_inf}'
        title.set_text(title_text)

        cbar.update_normal(im)
        cbar.set_label(label)

        draw_quiver(u, v, density_slider.val, scale_slider.val)
        fig.canvas.draw_idle()

    def on_field_change(label_text):
        nonlocal current_field
        current_field = label_text
        update()

    def on_vector_toggle(label_text):
        nonlocal show_vectors
        show_vectors = not show_vectors
        update()

    slider.on_changed(update)
    density_slider.on_changed(update)
    scale_slider.on_changed(update)
    radio.on_clicked(on_field_change)
    check.on_clicked(on_vector_toggle)

    plt.show()
    f.close()


if __name__ == '__main__':
    main()
