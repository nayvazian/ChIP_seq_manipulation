{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "id": "88e32cc9",
   "metadata": {},
   "outputs": [],
   "source": [
    "from pathlib import Path\n",
    "import pandas as pd\n",
    "\n",
    "EP300_PATH = Path(\"../src/EP300_ChIP-seq_peaks.bed\")\n",
    "CTCF_PATH = Path(\"../src/CTCF_ChIP-seq_peaks.bed\")\n",
    "K562_PATH = Path(\"../src/K562-Enhancers.bed\")\n",
    "\n",
    "def load_df(path):\n",
    "    return pd.read_csv(path, sep='\\t', header=None).sort_values(by=[0, 1])\n",
    "\n",
    "EP300_df = load_df(EP300_PATH)\n",
    "CTCF_df = load_df(CTCF_PATH)\n",
    "K562_df = load_df(K562_PATH)\n",
    "\n",
    "def find_intersection(dataset1, dataset2):\n",
    "    intersection = []\n",
    "    for range1 in dataset1:\n",
    "        for range2 in dataset2:\n",
    "            start = max(range1[0], range2[0])\n",
    "            end = min(range1[1], range2[1])\n",
    "            if start <= end:\n",
    "                intersection.append((start, end))\n",
    "    return intersection\n",
    "\n",
    "\n",
    "CHROM = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',\n",
    "       'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',\n",
    "       'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22',\n",
    "       'chrX']\n",
    "\n",
    "def intersect_by_chrom(dataset1, dataset2, dataset3=None):\n",
    "    chrom_intersection = []\n",
    "    for group in CHROM:\n",
    "        filtered_df1 = dataset1[dataset1[0] == group].iloc[:, [1, 2]]\n",
    "        filtered_df2 = dataset2[dataset2[0] == group].iloc[:, [1, 2]]\n",
    "        \n",
    "        overlap_list = find_intersection(filtered_df1.values.tolist(),\n",
    "                                          filtered_df2.values.tolist())\n",
    "\n",
    "        if dataset3 is not None:\n",
    "            filtered_df3 = dataset3[dataset3[0] == group].iloc[:, [1, 2]]\n",
    "            overlap_list = find_intersection(overlap_list,\n",
    "                                             filtered_df3.values.tolist())\n",
    "\n",
    "        chrom_intersection.append(overlap_list)\n",
    "\n",
    "    return chrom_intersection"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 3,
   "id": "f323b980",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "6884\n",
      "3027\n",
      "4119\n"
     ]
    },
    {
     "ename": "NameError",
     "evalue": "name 'CTCF_df' is not defined",
     "output_type": "error",
     "traceback": [
      "\u001b[1;31m---------------------------------------------------------------------------\u001b[0m",
      "\u001b[1;31mNameError\u001b[0m                                 Traceback (most recent call last)",
      "Cell \u001b[1;32mIn[3], line 10\u001b[0m\n\u001b[0;32m      8\u001b[0m total_length3 \u001b[38;5;241m=\u001b[39m \u001b[38;5;28msum\u001b[39m(\u001b[38;5;28mlen\u001b[39m(sublist) \u001b[38;5;28;01mfor\u001b[39;00m sublist \u001b[38;5;129;01min\u001b[39;00m output3)\n\u001b[0;32m      9\u001b[0m \u001b[38;5;28mprint\u001b[39m(total_length3)\n\u001b[1;32m---> 10\u001b[0m output4 \u001b[38;5;241m=\u001b[39m intersect_by_chrom(EP300_DF,\u001b[43mCTCF_df\u001b[49m,K562_DF)\n\u001b[0;32m     11\u001b[0m total_length4 \u001b[38;5;241m=\u001b[39m \u001b[38;5;28msum\u001b[39m(\u001b[38;5;28mlen\u001b[39m(sublist) \u001b[38;5;28;01mfor\u001b[39;00m sublist \u001b[38;5;129;01min\u001b[39;00m output4)\n\u001b[0;32m     12\u001b[0m \u001b[38;5;28mprint\u001b[39m(total_length4)\n",
      "\u001b[1;31mNameError\u001b[0m: name 'CTCF_df' is not defined"
     ]
    }
   ],
   "source": [
    "output = intersect_by_chrom(K562_df,EP300_df)\n",
    "total_length = sum(len(sublist) for sublist in output)\n",
    "print(total_length)\n",
    "output2 = intersect_by_chrom(K562_df,CTCF_df)\n",
    "total_length2 = sum(len(sublist) for sublist in output2)\n",
    "print(total_length2)\n",
    "output3 = intersect_by_chrom(EP300_df,CTCF_df)\n",
    "total_length3 = sum(len(sublist) for sublist in output3)\n",
    "print(total_length3)\n",
    "output4 = intersect_by_chrom(EP300_df,CTCF_df,K562_df)\n",
    "total_length4 = sum(len(sublist) for sublist in output4)\n",
    "print(total_length4)\n",
    "output4 = intersect_by_chrom(EP300_df,CTCF_df,K562_df)\n",
    "total_length4 = sum(len(sublist) for sublist in output4)\n",
    "print(total_length4)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "015af32e",
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.10.8"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
