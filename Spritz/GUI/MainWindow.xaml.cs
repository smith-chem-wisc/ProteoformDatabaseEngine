﻿using System;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Threading;
using System.Text.RegularExpressions;

namespace SpritzGUI
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        private readonly ObservableCollection<RNASeqFastqDataGrid> RnaSeqFastqCollection = new ObservableCollection<RNASeqFastqDataGrid>();
        private ObservableCollection<InRunTask> DynamicTasksObservableCollection = new ObservableCollection<InRunTask>();
        private readonly ObservableCollection<PreRunTask> StaticTasksObservableCollection = new ObservableCollection<PreRunTask>();
        private readonly ObservableCollection<SRADataGrid> SraCollection = new ObservableCollection<SRADataGrid>();
        private CancellationTokenSource TokenSource = new CancellationTokenSource();
        private EverythingRunnerEngine Everything;
        private Regex outputScrub = new Regex(@"(\[\d+m)");
        //private Task EverythingTask;

        public MainWindow()
        {
            InitializeComponent();
            DataGridRnaSeqFastq.DataContext = RnaSeqFastqCollection;
            workflowTreeView.DataContext = StaticTasksObservableCollection;
            LbxSRAs.ItemsSource = SraCollection;
            MessageBox.Show("Please have Docker Desktop installed and enable all shared drives.", "Setup", MessageBoxButton.OK, MessageBoxImage.Information);

            //var watch = new FileSystemWatcher();
            //watch.Path = Path.Combine(Environment.CurrentDirectory, "output");
            //watch.Filter = "test.txt";// Path.GetFileName(Everything.PathToWorkflow);
            //watch.NotifyFilter = NotifyFilters.LastWrite;
            //watch.Changed += new FileSystemEventHandler(OnWorkflowOutputChanged);
            //watch.EnableRaisingEvents = true;
        }

        protected override void OnClosed(EventArgs e)
        {
            // TODO: implement some way of killing EverythingTask

            // new process that kills docker container (if any)
            Process proc = new Process();
            proc.StartInfo.FileName = "Powershell.exe";
            proc.StartInfo.Arguments = "docker kill spritz";
            proc.StartInfo.CreateNoWindow = true;
            proc.StartInfo.UseShellExecute = false;
            proc.StartInfo.RedirectStandardError = true;
            proc.Start();

            if (proc != null && !proc.HasExited)
            {
                proc.WaitForExit();
            }

            base.OnClosed(e);
        }

        private void UpdateSRABox()
        {
            if (RnaSeqFastqCollection.Count > 0)
            {
                TbxSRA.IsEnabled = false;
                BtnAddSRA.IsEnabled = false;
                BtnClearSRA.IsEnabled = false;
            }
            else
            {
                TbxSRA.IsEnabled = true;
                BtnAddSRA.IsEnabled = true;
                BtnClearSRA.IsEnabled = true;
            }
        }

        private void Window_Drop(object sender, DragEventArgs e)
        {
            string[] files = (string[])e.Data.GetData(DataFormats.FileDrop);
            if (files != null)
            {
                foreach (var draggedFilePath in files)
                {
                    if (Directory.Exists(draggedFilePath))
                    {
                        foreach (string file in Directory.EnumerateFiles(draggedFilePath, "*.*", SearchOption.AllDirectories))
                        {
                            AddAFile(file);
                        }
                    }
                    else
                    {
                        AddAFile(draggedFilePath);
                    }
                    DataGridRnaSeqFastq.Items.Refresh();
                }
            }
            UpdateOutputFolderTextbox();
            UpdateSRABox();
        }

        private void Window_Loaded(object sender, RoutedEventArgs e)
        {
        }

        private void MenuItem_Wiki_Click(object sender, RoutedEventArgs e)
        {
            System.Diagnostics.Process.Start(@"https://github.com/smith-chem-wisc/Spritz/wiki");
        }

        private void MenuItem_Contact_Click(object sender, RoutedEventArgs e)
        {
            System.Diagnostics.Process.Start(@"https://github.com/smith-chem-wisc/Spritz");
        }

        private void RunWorkflowButton_Click(object sender, RoutedEventArgs e)
        {
            try
            {
                if (SraCollection.Count == 0 && RnaSeqFastqCollection.Count == 0)
                {
                    MessageBox.Show("You have not added any nucleic acid sequencing data (SRA accession or fastq files).", "Workflow", MessageBoxButton.OK, MessageBoxImage.Warning);
                    return;
                }

                if (StaticTasksObservableCollection.Count == 0)
                {
                    MessageBox.Show("You must add a workflow before a run.", "Run Workflows", MessageBoxButton.OK, MessageBoxImage.Information);
                    return;
                }
                else if (RnaSeqFastqCollection.Any() && GetPathToFastqs().CompareTo(OutputFolderTextBox.Text) != 0) // to be edited
                {
                    MessageBox.Show("FASTQ files do not exist in the user-defined analysis directory.", "Run Workflows", MessageBoxButton.OK, MessageBoxImage.Information);
                    return;
                }

                DynamicTasksObservableCollection = new ObservableCollection<InRunTask>();
                DynamicTasksObservableCollection.Add(new InRunTask("Workflow 1", StaticTasksObservableCollection.First().options));
                workflowTreeView.DataContext = DynamicTasksObservableCollection;
                
                Everything = new EverythingRunnerEngine(DynamicTasksObservableCollection.Select(b => new Tuple<string, Options>(b.DisplayName, b.options)).First(), OutputFolderTextBox.Text);

                WarningsTextBox.Document.Blocks.Clear();
                WarningsTextBox.AppendText($"Command executing: Powershell.exe {Everything.GenerateCommandsDry()}\n\n"); // keep for debugging
                WarningsTextBox.AppendText($"Saving output to {Everything.PathToWorkflow}. Please monitor it there...\n\n");

                Everything.WriteConfig(StaticTasksObservableCollection.First().options);
                var t = new Task(RunEverythingRunner);
                t.Start();
                t.ContinueWith(DisplayAnyErrors);

                // update gui
                RunWorkflowButton.IsEnabled = false;
                ClearTasksButton.IsEnabled = true;
                BtnWorkFlow.IsEnabled = false;
                ResetTasksButton.IsEnabled = true;
            }
            catch (TaskCanceledException)
            {
                // Ignore error
            }
        }

        private void RunEverythingRunner()
        {
            Process proc = new Process();
            proc.StartInfo.FileName = "Powershell.exe";
            proc.StartInfo.Arguments = Everything.GenerateCommandsDry();
            proc.StartInfo.UseShellExecute = false;
            proc.StartInfo.RedirectStandardOutput = true;
            proc.StartInfo.RedirectStandardError = true;
            proc.StartInfo.CreateNoWindow = true;
            proc.OutputDataReceived += new DataReceivedEventHandler(OutputHandler);
            proc.ErrorDataReceived += new DataReceivedEventHandler(OutputHandler);
            proc.Start();
            proc.BeginOutputReadLine();
            proc.BeginErrorReadLine();
            proc.WaitForExit();
        }

        private void OutputHandler(object source, DataReceivedEventArgs e)
        {
            Dispatcher.Invoke(() => 
            {
                string output = outputScrub.Replace(e.Data, "");
                WarningsTextBox.AppendText(output + Environment.NewLine);
                using (StreamWriter sw = File.Exists(Everything.PathToWorkflow) ? File.AppendText(Everything.PathToWorkflow) : File.CreateText(Everything.PathToWorkflow))
                {
                    sw.WriteLine(output);
                }
            });
        }

        private void DisplayAnyErrors(Task obj)
        {
            Dispatcher.Invoke(() => WarningsTextBox.AppendText("Done!" + Environment.NewLine));
            Dispatcher.Invoke(() => MessageBox.Show("Finished! Workflow summary is located in " 
                + StaticTasksObservableCollection.First().options.AnalysisDirectory, "Spritz Workflow", 
                MessageBoxButton.OK, MessageBoxImage.Information));
        }

        private void BtnAddRnaSeqFastq_Click(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.OpenFileDialog openPicker = new Microsoft.Win32.OpenFileDialog()
            {
                Filter = "FASTQ Files|*.fastq",
                FilterIndex = 1,
                RestoreDirectory = true,
                Multiselect = true
            };
            if (openPicker.ShowDialog() == true)
            {
                foreach (var filepath in openPicker.FileNames)
                {
                    AddAFile(filepath);
                }
            }
            DataGridRnaSeqFastq.Items.Refresh();
            UpdateSRABox();
        }

        private void BtnClearRnaSeqFastq_Click(object sender, RoutedEventArgs e)
        {
            RnaSeqFastqCollection.Clear();
            UpdateOutputFolderTextbox();
            UpdateSRABox();
        }

        //private void LoadTaskButton_Click(object sender, RoutedEventArgs e)
        //{
        //    Microsoft.Win32.OpenFileDialog openPicker = new Microsoft.Win32.OpenFileDialog()
        //    {
        //        Filter = "TOML files(*.toml)|*.toml",
        //        FilterIndex = 1,
        //        RestoreDirectory = true,
        //        Multiselect = true
        //    };
        //    if (openPicker.ShowDialog() == true)
        //    {
        //        foreach (var tomlFromSelected in openPicker.FileNames)
        //        {
        //            AddAFile(tomlFromSelected);
        //        }
        //    }
        //    UpdateTaskGuiStuff();
        //}

        private void ClearTasksButton_Click(object sender, RoutedEventArgs e)
        {
            StaticTasksObservableCollection.Clear();
            workflowTreeView.DataContext = StaticTasksObservableCollection;
            WarningsTextBox.Document.Blocks.Clear();
            UpdateTaskGuiStuff();
        }

        private void ResetTasksButton_Click(object sender, RoutedEventArgs e)
        {
            RunWorkflowButton.IsEnabled = true;
            ClearTasksButton.IsEnabled = true;
            BtnWorkFlow.IsEnabled = false;
            ResetTasksButton.IsEnabled = false;

            DynamicTasksObservableCollection.Clear();
            workflowTreeView.DataContext = StaticTasksObservableCollection;
        }

        //private void AddNewRnaSeqFastq(object sender, StringListEventArgs e)
        //{
        //    if (!Dispatcher.CheckAccess())
        //    {
        //        Dispatcher.BeginInvoke(new Action(() => AddNewRnaSeqFastq(sender, e)));
        //    }
        //    else
        //    {
        //        foreach (var uu in RnaSeqFastqCollection)
        //        {
        //            uu.Use = false;
        //        }
        //        foreach (var newRnaSeqFastqData in e.StringList)
        //        {
        //            RnaSeqFastqCollection.Add(new RNASeqFastqDataGrid(newRnaSeqFastqData));
        //        }
        //        UpdateOutputFolderTextbox();
        //    }
        //}

        private void BtnAddSRA_Click(object sender, RoutedEventArgs e)
        {
            if (TbxSRA.Text.Contains("SR") || TbxSRA.Text.Contains("ER"))
            {
                if (SraCollection.Any(s => s.Name == TbxSRA.Text.Trim()))
                {
                    MessageBox.Show("That SRA has already been added. Please choose a new SRA accession.", "Workflow", MessageBoxButton.OK, MessageBoxImage.Information);
                }
                else
                {
                    SRADataGrid sraDataGrid = new SRADataGrid(TbxSRA.Text.Trim());
                    SraCollection.Add(sraDataGrid);
                }
            }
            else if (MessageBox.Show("SRA accessions are expected to start with \"SR\" or \"ER\", such as SRX254398 or ERR315327. View the GEO SRA website?", "Workflow", MessageBoxButton.YesNo, MessageBoxImage.Question, MessageBoxResult.No) == MessageBoxResult.Yes)
            {
                System.Diagnostics.Process.Start("https://www.ncbi.nlm.nih.gov/sra");
            }
        }

        private void BtnClearSRA_Click(object sender, RoutedEventArgs e)
        {
            SraCollection.Clear();
            BtnAddSRA.IsEnabled = true;
        }

        private void BtnWorkFlow_Click(object sender, RoutedEventArgs e)
        {
            if (SraCollection.Count == 0 && RnaSeqFastqCollection.Count == 0)
            {
                if (MessageBox.Show("You have not added any nucleic acid sequencing data (SRA accession or fastq files). Would you like to continue to make a protein database from the reference gene model?", "Workflow", MessageBoxButton.YesNo, MessageBoxImage.Question, MessageBoxResult.No) == MessageBoxResult.No)
                {
                    return;
                }
            }

            try
            {
                var dialog = new WorkFlowWindow(OutputFolderTextBox.Text == "" ? new Options().AnalysisDirectory : OutputFolderTextBox.Text);
                if (dialog.ShowDialog() == true)
                {
                    AddTaskToCollection(dialog.Options);
                    UpdateTaskGuiStuff();
                    UpdateOutputFolderTextbox();
                }
            }
            catch (InvalidOperationException)
            {
                // does not open workflow window until all fastq files are added, if any
            }
        }

        //private void BtnSaveRnaSeqFastqSet_Click(object sender, RoutedEventArgs e)
        //{
        //    try
        //    {
        //        WriteExperDesignToTsv(OutputFolderTextBox.Text);
        //    }
        //    catch (Exception ex)
        //    {
        //        MessageBox.Show("Could not save experimental design!\n\n" + ex.Message, "Experimental Design", MessageBoxButton.OK, MessageBoxImage.Warning);
        //        return;
        //    }
        //}

        private void UpdateTaskGuiStuff()
        {
            if (StaticTasksObservableCollection.Count == 0)
            {
                RunWorkflowButton.IsEnabled = false;
                ClearTasksButton.IsEnabled = false;
                BtnWorkFlow.IsEnabled = true;
                ResetTasksButton.IsEnabled = false;
            }
            else
            {
                RunWorkflowButton.IsEnabled = true;
                ClearTasksButton.IsEnabled = true;
                BtnWorkFlow.IsEnabled = false;
                ResetTasksButton.IsEnabled = false;
            }
        }

        private void AddTaskToCollection(Options ye)
        {
            PreRunTask te = new PreRunTask(ye);
            StaticTasksObservableCollection.Add(te);
            StaticTasksObservableCollection.Last().DisplayName = "Task" + (StaticTasksObservableCollection.IndexOf(te) + 1);
        }

        private string GetPathToFastqs()
        {
            var MatchingChars =
                    from len in Enumerable.Range(0, RnaSeqFastqCollection.Select(b => b.FilePath).Min(s => s.Length)).Reverse()
                    let possibleMatch = RnaSeqFastqCollection.Select(b => b.FilePath).First().Substring(0, len)
                    where RnaSeqFastqCollection.Select(b => b.FilePath).All(f => f.StartsWith(possibleMatch, StringComparison.Ordinal))
                    select possibleMatch;

            return Path.Combine(Path.GetDirectoryName(MatchingChars.First()));
        }

        //private string GetPathToFastqDirectory(string path)
        //{
        //    var filePath = path.Split('\\');
        //    var newPath = string.Join("\\", filePath.Take(filePath.Length - 1));
        //    return newPath;
        //}

        //private void UpdateOutputFolderTextbox(string filePath = null)
        //{
        //    // if new files have different path than current text in output, then throw error
        //    if (StaticTasksObservableCollection.Count > 0)
        //    {
        //        OutputFolderTextBox.Text = StaticTasksObservableCollection.First().options.AnalysisDirectory;
        //    }
        //    else if (RnaSeqFastqCollection.Any())
        //    {
        //        if (filePath != null && OutputFolderTextBox.Text != "" && GetPathToFastqDirectory(filePath).CompareTo(OutputFolderTextBox.Text) != 0)
        //        {
        //            throw new InvalidOperationException();
        //        }
        //        OutputFolderTextBox.Text = GetPathToFastqs();
        //    }
        //    else
        //    {
        //        OutputFolderTextBox.Clear();
        //    }
        //}

        private void UpdateOutputFolderTextbox()
        {
            if (StaticTasksObservableCollection.Count > 0)
            {
                OutputFolderTextBox.Text = StaticTasksObservableCollection.First().options.AnalysisDirectory;
            }
            else if (RnaSeqFastqCollection.Any())
            {
                OutputFolderTextBox.Text = GetPathToFastqs();
            }
            else
            {
                OutputFolderTextBox.Clear();
            }
        }

        private void AddAFile(string filepath)
        {
            if (SraCollection.Count == 0)
            {
                var theExtension = Path.GetExtension(filepath).ToLowerInvariant();
                theExtension = theExtension == ".gz" ? Path.GetExtension(Path.GetFileNameWithoutExtension(filepath)).ToLowerInvariant() : theExtension;
                switch (theExtension)
                {
                    case ".fastq":
                        if (Path.GetFileName(filepath).Contains("_1") || Path.GetFileName(filepath).Contains("_2"))
                        {
                            RNASeqFastqDataGrid rnaSeqFastq = new RNASeqFastqDataGrid(filepath);
                            RnaSeqFastqCollection.Add(rnaSeqFastq);
                            UpdateOutputFolderTextbox();
                            break;
                        }
                        else
                        {
                            MessageBox.Show("FASTQ files must have *_1.fastq and *_2.fastq extensions.", "Run Workflows", MessageBoxButton.OK, MessageBoxImage.Information);
                            return;
                        }
                }
            }
            else
            {
                MessageBox.Show("User already added SRA number. Please only choose one input: 1) SRA accession 2) FASTQ files.", "Run Workflows", MessageBoxButton.OK, MessageBoxImage.Information);
                return;
            }
        }

        private void workflowTreeView_MouseDoubleClick(object sender, MouseButtonEventArgs e)
        {
            var a = sender as TreeView;
            if (a.SelectedItem is PreRunTask preRunTask)
            {
                var workflowDialog = new WorkFlowWindow(preRunTask.options);
                workflowDialog.ShowDialog();
                workflowTreeView.Items.Refresh();
            }
        }

        private void WarningsTextBox_TextChanged(object sender, TextChangedEventArgs e)
        {
            WarningsTextBox.ScrollToEnd();
        }
    }
}